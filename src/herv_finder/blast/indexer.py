"""A general-purposed multi-processed BLAST index creator"""
import copy
import gc
import logging
import multiprocessing
import queue
import signal
import time
from collections import defaultdict
from typing import Tuple, Dict, Optional, Set, Union

import tqdm

from herv_finder.utils import compressed_pickle
from herv_finder.utils import in_memory_fasta
from herv_finder.utils import parallel_helper

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()

blast_index_data_type = Set[Tuple[str, bool, int]]

blast_index_type = Dict[bytes, blast_index_data_type]
"""
The type of blast index. Will be a Dictionary with:

index_item -> [chromosome_name, strand, offset]

strand is True for default and False for reversed.
"""

DEFAULT_INDEX_LEN = 4
"""Index length of nucleotides"""

DEFAULT_CHUNK_LEN = 2000000
"""Default chunk length"""


def _blast_index_merge(d1: blast_index_type, d2: blast_index_type) -> blast_index_type:
    dd1 = defaultdict(lambda: set())
    dd1.update(d1)
    for k, v in d2.items():
        dd1[k].update(v)
    return copy.deepcopy(dict(dd1))


def _is_low_complexity(fasta_bytes: bytes) -> bool:
    """
    Whether the sequence is of low complexity.

    TODO
    """
    # return b'N' in fasta_bytes or b'a' in fasta_bytes
    _ = fasta_bytes
    return False


class _IndexWorkerProcess(multiprocessing.Process):
    def __init__(self,
                 chromosome_name: str,
                 strand: bool,
                 fasta_bytes: bytes,
                 index_len: int,
                 start_pos: int,
                 tmp_dir: multiprocessing.Queue
                 ):
        """
        The worker of index creator. Will create index for both forward and reverse strand.

        :param chromosome_name: Name of chromosome.
        :param strand: Whether this is + strand.
        :param fasta_bytes: Bytes needed to be indexed.
        :param index_len: The length (k of the k-mer) of the index.
        :param start_pos: 0-based offset of the fasta_bytes.
        :param tmp_dir: Temporary directory.
        """
        super().__init__()
        self.start_pos = start_pos
        self.index_len = index_len
        self.fasta_bytes = fasta_bytes
        self.strand = strand
        self.chromosome_name = chromosome_name
        self.tmp_dir = tmp_dir
        self.logger_handler = logging.getLogger()
        self.logger_handler.debug(self.__repr__())
        self.tmp_dict = defaultdict(lambda: set())


    def run(self):
        self.logger_handler.debug(f"Worker PID: {self.pid} started")
        for bytes_window_start in range(len(self.fasta_bytes) - self.index_len + 1):
            # print(self.start_pos+bytes_window_start)
            bytes_window = self.fasta_bytes[bytes_window_start:bytes_window_start + self.index_len]
            if not self.strand:
                bytes_window = in_memory_fasta.get_reversed_complementary(bytes_window)
            if _is_low_complexity(bytes_window):
                continue
            # if not bytes_window in self.tmp_dict:
            #     self.tmp_dict[bytes_window] = set()
            self.tmp_dict[bytes_window].add(
                (self.chromosome_name, self.strand, self.start_pos + bytes_window_start)
            )
        self.tmp_dir.put(dict(self.tmp_dict))

    def __repr__(self):
        try:
            return f"Worker {self.chromosome_name} " + \
                   f"[{self.strand}] {self.start_pos} -> {self.start_pos + len(self.fasta_bytes)}"
        except AttributeError:
            return "Worker under construction"


class _MergeWorkerProcess(multiprocessing.Process):
    def __init__(self,
                 input_queue: multiprocessing.Queue,
                 output_list: list,
                 supervising_pool: parallel_helper.MultiProcessingJobQueue
                 ):
        super().__init__()
        self.logger_handler = logging.getLogger("Merger")
        self.input_queue = input_queue
        self.output_list = output_list
        self.multiprocessing = supervising_pool
        self.merged_dict = {}
        self.is_terminated = False
        signal.signal(signal.SIGTERM, self.sigterm)

    def run(self):
        time.sleep(1)
        self.logger_handler.info("Started")
        while not self.is_terminated:
            try:
                get_item = self.input_queue.get(block=False, timeout=0.1)
            except queue.Empty:
                continue
            except TimeoutError:
                continue
            self.logger_handler.info("Got finished process! Start merging...")
            self.merged_dict = _blast_index_merge(get_item, self.merged_dict)
            del get_item
        self.logger_handler.info("FIN recv, transmitting...")
        self.output_list.append(self.merged_dict)
        self.logger_handler.info("FIN")

    def sigterm(self, *args, **kwargs):
        self.is_terminated = True

class BlastIndex:
    def __init__(self, pool_len: Optional[int] = multiprocessing.cpu_count()):
        logger_handler.info("BlastIndex initializing...")
        self._indicies = []
        """A list of unmerged indicies"""

        self._pool_len = pool_len
        """Internal thread number"""

        self._fasta_obj = None
        """Internal Fasta object"""

        self.index_len = None

        self.total_index_count = 0
        """Number of entries inside the indicies. Can be duplicated."""

    def create_index(self, fasta_filename: str, index_len: Optional[int] = DEFAULT_INDEX_LEN):
        """
        Wrapper function to (re)-create the index.

        :param fasta_filename: Filename of input FASTA
        :param index_len: The length (k of the k-mer) of the index.
        """
        logger_handler.info(f"Creating index from {fasta_filename}...")
        self._fasta_obj = in_memory_fasta.Fasta(fasta_filename)

        self.index_len = index_len
        self._indicies = []

        pool = parallel_helper.MultiProcessingJobQueue(
            pool_name="Creating chunk-level Index",
            pool_size=self._pool_len,
            with_tqdm=True,
        )
        """The pool of process"""

        sync_manager = multiprocessing.Manager()
        public_index = sync_manager.list()
        public_queue = queue.Queue(maxsize=-1)
        # Do not use multiprocessing.Queue.
        new_worker_merger = _MergeWorkerProcess(public_queue, public_index, pool)
        new_worker_merger.start()
        for chromosome_name in self._fasta_obj.chromosomes:
            """Submit index creation job for a particular chromosome."""
            # The forward strand
            fasta_bytes = self._fasta_obj.get(chromosome_name)
            """The chromosome in bytes"""

            total_len = len(fasta_bytes)
            """Total length of this chromosome"""

            chunk_len = DEFAULT_CHUNK_LEN
            """Length of each chunk"""

            # Submit the jobs
            for i in range(total_len // chunk_len + 1):
                input_fasta_bytes = fasta_bytes[
                                    chunk_len * i:chunk_len * (i + 1) + self.index_len - 1
                                    ]
                new_worker = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                                 strand=True,
                                                 fasta_bytes=input_fasta_bytes,
                                                 index_len=self.index_len,
                                                 start_pos=chunk_len * i,
                                                 tmp_dir=public_queue)
                pool.append(new_worker)
                new_worker_rc = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                                    strand=False,
                                                    fasta_bytes=input_fasta_bytes,
                                                    index_len=self.index_len,
                                                    start_pos=chunk_len * i,
                                                    tmp_dir=public_queue)
                pool.append(new_worker_rc)

        pool.start()
        pool.join()
        new_worker_merger.terminate()
        new_worker_merger.join()


        self._indicies = public_index
        del public_index
        new_worker_merger.close()
        gc.collect()
        self._update_total_index_count()

    def show_index(self) -> str:
        rets = ""
        for index in self._indicies:
            for k, v in index.items():
                rets += f"{k}: {v}\n"
        return rets

    def write_index(self, dest_filename: str):
        """Write index to a file."""
        logger_handler.info(f"Writing index to {dest_filename}...")
        compressed_pickle.dump(self._indicies, dest_filename)
        logger_handler.info(f"Writing index to {dest_filename} FIN")

    def read_index(self, from_filename: str):
        """Read index from a file."""
        logger_handler.info(f"Reading index from {from_filename}...")
        self._indicies = compressed_pickle.load_with_tqdm(from_filename)
        self.index_len = len(list(self._indicies[0].keys())[0])
        logger_handler.info(f"Reading index from {from_filename} FIN")

    def get_stats(self):
        """Print statistics of this index"""
        print(f"Index length: {self.index_len}")
        self._update_total_index_count()
        value_len_dir = defaultdict(lambda: 0)
        value_len_sum = 0
        for index in self._indicies:
            for v in index.values():
                v_len = len(v)
                value_len_dir[v_len] += 1
            value_len_sum += sum(value_len_dir.values())
        print(
            f"Value len distribution: {value_len_dir}, sum {value_len_sum}, mean {value_len_sum / self.total_index_count}")

    def valid_index(self):
        """Validity the index. Should NOT be called by outer process."""
        if self._fasta_obj is None:
            return True

        with tqdm.tqdm(desc="Validating index", total=self.total_index_count) as pbar:
            for index in self._indicies:
                for k, v in index.items():
                    for location in v:
                        assert self._fasta_obj.get(
                            chromosome_name=location[0],
                            start=location[2],
                            end=location[2] + self.index_len,
                            strand=location[1]) == k
                        pbar.update(1)

    def get_pos(self, fasta_bytes: bytes) -> blast_index_data_type:
        rets = set()
        for index in self._indicies:
            try:
                rets.update(index[fasta_bytes])
            except KeyError:
                pass
        return rets

    def _update_total_index_count(self):
        self.total_index_count = 0
        for index in self._indicies:
            for v in index.values():
                self.total_index_count += len(v)


if __name__ == '__main__':
    bi = BlastIndex()
    bi.create_index('test/e_coli.fasta')
    # print(bi.show_index())
    bi.valid_index()
    bi.get_stats()
    # bi.write_index('test/test.pkl.xz')
    # bj = BlastIndex()
    # bj.read_index('test/test.pkl.xz')
