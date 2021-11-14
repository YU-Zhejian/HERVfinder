"""A general-purposed multi-processed BLAST index creator"""
import copy
import gc
import glob
import logging
import multiprocessing
import os
import pickle
import queue
import signal
import tempfile
import threading
import time
from collections import defaultdict
from typing import Optional, List

import tqdm

from herv_finder.blast import blast_index_data_type, blast_index_type, DEFAULT_INDEX_LEN, DEFAULT_CHUNK_LEN
from herv_finder.utils import compressed_pickle
from herv_finder.utils import in_memory_fasta
from herv_finder.utils import parallel_helper

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()


def _blast_index_merge(d1: blast_index_type, d2: blast_index_type) -> blast_index_type:
    dd1 = defaultdict(lambda: [])
    dd1.update(d1)
    for k, v in d2.items():
        dd1[k].extend(v)
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
                 fasta_bytes: bytes,
                 index_len: int,
                 start_pos: int,
                 output_queue: multiprocessing.Queue
                 ):
        """
        The worker of index creator. Will create index for both forward and reverse strand.

        :param chromosome_name: Name of chromosome.
        :param fasta_bytes: Bytes needed to be indexed.
        :param index_len: The length (k of the k-mer) of the index.
        :param start_pos: 0-based offset of the fasta_bytes.
        :param output_queue: Temporary directory.
        """
        super().__init__()
        self.start_pos = start_pos
        self.index_len = index_len
        self.fasta_bytes = fasta_bytes
        self.chromosome_name = chromosome_name
        self.output_queue = output_queue
        self.logger_handler = logging.getLogger("Worker")
        self.logger_handler.debug(self.__repr__())
        self.tmp_dict = None

    def run(self):
        self.tmp_dict = defaultdict(lambda: [])
        self.logger_handler.debug(f"Worker PID: {self.pid} started")
        for bytes_window_start in range(len(self.fasta_bytes) - self.index_len + 1):
            # print(self.start_pos+bytes_window_start)
            bytes_window = self.fasta_bytes[bytes_window_start:bytes_window_start + self.index_len]
            if not _is_low_complexity(bytes_window):
                self.tmp_dict[bytes_window].append(
                    (self.chromosome_name, True, self.start_pos + bytes_window_start)
                )
            bytes_window_rc = in_memory_fasta.get_reversed_complementary(bytes_window)
            if not _is_low_complexity(bytes_window_rc):
                self.tmp_dict[bytes_window_rc].append(
                    (self.chromosome_name, False, self.start_pos + bytes_window_start)
                )
        self.logger_handler.debug(f"Worker PID: {self.pid} start transmitting data")
        self.output_queue.put(dict(self.tmp_dict), block=False)
        self.logger_handler.debug(f"Worker PID: {self.pid} FIN")

    def __repr__(self):
        try:
            return f"Worker {self.chromosome_name} " + \
                   f"{self.start_pos} -> {self.start_pos + len(self.fasta_bytes)}"
        except AttributeError:
            return "Worker under construction"


class _MergeWorkerProcess(multiprocessing.Process):
    """Merge the chunk-level indices and split it by two-mers."""

    def __init__(self,
                 input_queue: multiprocessing.Queue,
                 output_queue: multiprocessing.Queue,
                 per_chunk_index_finish_flag: threading.Event
                 ):
        super().__init__()
        self.logger_handler = logging.getLogger("Merger")
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.per_chunk_index_finish_flag = per_chunk_index_finish_flag

    def run(self):
        time.sleep(1)
        self.logger_handler.debug("Started")
        total_items_get = 0
        tmp_dir = tempfile.mkdtemp()
        while not self.per_chunk_index_finish_flag.is_set() or not self.input_queue.empty():
            try:
                get_item = self.input_queue.get(block=False, timeout=0.1)
            except queue.Empty:
                continue
            except TimeoutError:
                continue
            self.logger_handler.debug("Got finished process! Start merging...")
            total_items_get += 1
            for two_mer in in_memory_fasta.get_all_kmers(2, bases=b'AGCTNX'):
                two_mer_str = str(two_mer,encoding='UTF-8')
                index_started_by_two_mer = {}
                for k, v in get_item.items():
                    if k.startswith(two_mer):
                        index_started_by_two_mer[k] = v
                if len(index_started_by_two_mer) > 0:
                    output_filename = f"{two_mer_str}-{total_items_get}.pkl"
                    pickle.dump(index_started_by_two_mer, open(os.path.join(tmp_dir,output_filename), "wb"))
                    del index_started_by_two_mer
            del get_item
            gc.collect()
        self.logger_handler.info("Re-splitting the index to 2-mers...")
        for two_mer in tqdm.tqdm(desc="Merging...",
                                 iterable=list(in_memory_fasta.get_all_kmers(2, bases=b'AGCTNX'))
                                 ):
            two_mer_str = str(two_mer, encoding='UTF-8')
            index_started_by_two_mer = {}
            for filename in glob.glob(os.path.join(tmp_dir,f"{two_mer_str}-*.pkl")):
                self.logger_handler.debug(f"Joining chunk-level index from {filename}")
                get_item = pickle.load(open(filename,"rb"))
                index_started_by_two_mer=_blast_index_merge(index_started_by_two_mer,get_item)
                del get_item
            if len(index_started_by_two_mer) > 0:
                self.output_queue.put(
                    (two_mer, index_started_by_two_mer)
                )
                # del index_started_by_two_mer
        self.output_queue.close()
        self.logger_handler.debug("FIN")


class BlastIndex:
    """The blast index creator and utilities"""
    def __init__(self,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 index_len: Optional[int] = DEFAULT_INDEX_LEN,
                 chunk_len: Optional[int] = DEFAULT_CHUNK_LEN
                 ):
        logger_handler.info("BlastIndex initializing...")
        self._indices = {}
        """
        A dict of re-split indices. Format:
        
        Dict[leading-bytes: blast_index_type]
        """

        self.pool_len = pool_len
        """Internal thread number"""

        self._fasta_obj = None
        """Internal Fasta object"""

        self.index_len = index_len
        self.chunk_len = chunk_len

        self.total_index_count = 0
        """Number of entries inside the indicies. Can be duplicated."""

    def attach_fasta(self, fasta_filename: str):
        self._fasta_obj = in_memory_fasta.Fasta(fasta_filename)

    def detach_fasta(self):
        if self._fasta_obj is not None:
            del self._fasta_obj
            self._fasta_obj=None

    def create_index_singlethread(self):
        tmp_dict = defaultdict(lambda: [])

        self._indices = {}
        with tqdm.tqdm(total=2 * self._fasta_obj.total_length) as pbar:
            for chromosome_name in self._fasta_obj.chromosomes:
                """Submit index creation job for a particular chromosome."""
                # The forward strand
                fasta_bytes = self._fasta_obj.get(chromosome_name)

                for bytes_window_start in range(len(fasta_bytes) - self.index_len + 1):
                    # print(self.start_pos+bytes_window_start)
                    bytes_window = fasta_bytes[bytes_window_start:bytes_window_start + self.index_len]

                    if not _is_low_complexity(bytes_window):
                        tmp_dict[bytes_window].append(
                            (chromosome_name, True, bytes_window_start)
                        )
                    pbar.update(1)
                    bytes_window_rc = in_memory_fasta.get_reversed_complementary(bytes_window)
                    if not _is_low_complexity(bytes_window_rc):
                        tmp_dict[bytes_window_rc].append(
                            (chromosome_name, False, bytes_window_start)
                        )
                    pbar.update(1)
        self._indices = [dict(tmp_dict)] # FIXME
        self._update_total_index_count()

    def create_index(self, dest_basename: str):
        """
        Wrapper function to (re)-create the index.

        :param fasta_filename: Filename of input FASTA
        :param index_len: The length (k of the k-mer) of the index.
        """

        self._indices = {}

        pool = parallel_helper.MultiProcessingJobQueue(
            pool_name="Creating chunk-level Index",
            pool_size=self.pool_len,
            with_tqdm=True,
        )
        """The pool of process"""

        sync_manager = multiprocessing.Manager()
        per_prefix_index_queue = sync_manager.Queue(maxsize=-1)
        per_chunk_index_queue = sync_manager.Queue(maxsize=-1)
        per_chunk_index_finish_flag = sync_manager.Event()
        per_chunk_index_finish_flag.clear()
        # Do not use multiprocessing.Queue.
        worker_merger = _MergeWorkerProcess(
            input_queue=per_chunk_index_queue,
            output_queue=per_prefix_index_queue,
            per_chunk_index_finish_flag=per_chunk_index_finish_flag
        )
        worker_merger.start()
        for chromosome_name in self._fasta_obj.chromosomes:
            # Submit index creation job for a particular chromosome
            fasta_bytes = self._fasta_obj.get(chromosome_name)
            """The chromosome in bytes"""

            total_len = len(fasta_bytes)
            """Total length of this chromosome"""

            # Submit the jobs
            for i in range(total_len // self.chunk_len + 1):
                start_pos = self.chunk_len * i
                input_fasta_bytes = fasta_bytes[
                                    start_pos:
                                    start_pos + self.chunk_len + self.index_len - 1
                                    ]
                new_worker = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                                 fasta_bytes=input_fasta_bytes,
                                                 index_len=self.index_len,
                                                 start_pos=start_pos,
                                                 output_queue=per_chunk_index_queue)
                pool.append(new_worker)

        pool.start()
        pool.join()
        logger_handler.info(f"Chunk-level indexing had finished, re-merging in prefix-level...")
        per_chunk_index_finish_flag.set() # Notify the merger that all chunks have finished
        while worker_merger.exitcode is None:
            try:
                get_item = per_prefix_index_queue.get(block=False, timeout=0.1)
            except queue.Empty:
                continue
            except TimeoutError:
                continue
            kmer_str = str(get_item[0], encoding='UTF-8')
            output_filename = f"{dest_basename}_{kmer_str}.pkl.xz"
            logger_handler.debug(f"Writing index to {output_filename}...")
            compressed_pickle.dump(get_item[1], output_filename)
            del get_item
        print(worker_merger.exitcode)
        worker_merger.join()
        worker_merger.close()
        gc.collect()
        logger_handler.info(f"Index FIN")

    def show_index(self) -> str:
        rets = ""
        for two_mer, index in self._indices.items():
            rets += str(two_mer, encoding='UTF-8') + '\n'
            for k, v in index.items():
                rets += f"{k}: {v}\n"
            rets += '\n'
        return rets

    def read_index(self, from_basename: str, two_mers:List[bytes] = in_memory_fasta.get_all_kmers(2)):
        """Read index from a file."""
        logger_handler.info(f"Reading index from {from_basename}...")
        self._indices ={}
        for two_mer in two_mers:
            kmer_str = str(two_mer, encoding='UTF-8')
            load_filename = f"{from_basename}_{kmer_str}.pkl.xz"
            if os.path.exists(load_filename):
                logger_handler.info(f"Reading index from {load_filename}...")
                self._indices[two_mer] = compressed_pickle.load_with_tqdm(load_filename)
        try:
            self.index_len = len(list(list(self._indices.values())[0].keys())[0])
        except:
            self.index_len = None
        logger_handler.info(f"Reading index from {from_basename} FIN")

    def get_stats(self):
        """Print statistics of this index"""
        print(f"Index length: {self.index_len}")
        self._update_total_index_count()
        if self.total_index_count == 0:
            return
        value_len_dir = defaultdict(lambda: 0)
        value_len_sum = 0
        for index in self._indices.values():
            for v in index.values():
                v_len = len(v)
                value_len_dir[v_len] += 1
                value_len_sum += v_len
        print(
            f"Value len distribution: sum {value_len_sum}, mean {value_len_sum / self.total_index_count}",
            str(dict(value_len_dir))
        )

    def valid_index(self):
        """Validity the index. Should NOT be called by outer process."""
        if self._fasta_obj is None:
            return

        with tqdm.tqdm(desc="Validating index", total=self.total_index_count) as pbar:
            for index in self._indices.values():
                for k, v in index.items():
                    for location in v:
                        assert self._fasta_obj.get(
                            chromosome_name=location[0],
                            start=location[2],
                            end=location[2] + self.index_len,
                            strand=location[1]) == k
                        pbar.update(1)

    def get_pos_singlethread(self, fasta_bytes: bytes) -> blast_index_data_type:
        rets = []
        for index in self._indices.values():
            try:
                rets.extend(index[fasta_bytes])
            except KeyError:
                pass
        return rets

    def get_pos(self, fasta_bytes: bytes):
        pass

    @property
    def keys(self):
        return self._indices.keys()

    @property
    def keys_full_length(self):
        return [v.keys() for v in self._indices.values()]

    def _update_total_index_count(self):
        self.total_index_count = 0
        for index in tqdm.tqdm(desc="Updating total index count",
                               iterable=self._indices.values()):
            for v in index.values():
                self.total_index_count += len(v)


def _test_on_tiny():
    bi = BlastIndex(index_len=4, chunk_len=20)
    bi.attach_fasta('test/tiny.fasta')
    bi.create_index('test/test')
    bi.read_index('test/test')
    print(bi.show_index())
    bi.valid_index()
    bi.detach_fasta()
    bi.get_stats()
    del bi
    bj = BlastIndex()
    bj.read_index('test/test')
    print(bj.show_index())
    bj.valid_index()
    bj.get_stats()


if __name__ == '__main__':
    _test_on_tiny()
