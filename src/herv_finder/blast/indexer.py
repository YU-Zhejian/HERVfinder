"""A general-purposed multi-processed BLAST index creator"""
import copy
import glob
import logging
import multiprocessing
import os
import pickle
import tempfile
import threading
from collections import defaultdict
from typing import Tuple, Dict, Optional, Set

import concurrent.futures
import tqdm

from herv_finder.utils import compressed_pickle
from herv_finder.utils import parallel_helper
from herv_finder.utils import in_memory_fasta

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()

blast_index_type = Dict[bytes, Set[Tuple[str, bool, int]]]
"""
The type of blast index. Will be a Dictionary with:

index_item -> [chromosome_name, strand, offset]

strand is True for default and False for reversed.
"""

DEFAULT_INDEX_LEN = 11
"""Index length of nucleotides"""

DEFAULT_CHUNK_LEN = 1000000
"""Default chunk length"""


def _blast_index_merge(d1: blast_index_type, d2: blast_index_type) -> blast_index_type:
    dd1 = defaultdict(lambda: set())
    dd1.update(d1)
    for k, v in d2.items():
        dd1[k].update(v)
    return dict(dd1)


def _is_low_complexity(fasta_bytes: bytes) -> bool:
    """
    Whether the sequence is of low complexity.

    TODO
    """
    # return b'N' in fasta_bytes or b'a' in fasta_bytes
    _=fasta_bytes
    return False


class _IndexWorkerProcess(multiprocessing.Process):
    def __init__(self,
                 chromosome_name: str,
                 strand: bool,
                 fasta_bytes: bytes,
                 index_len: int,
                 start_pos: int,
                 tmp_dir: str
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
        if not self.strand:
            self.fasta_bytes= in_memory_fasta.get_reversed_complementary(self.fasta_bytes)
        self.chromosome_name = chromosome_name
        self.tmp_dir = tmp_dir
        self.logger_handler = logging.getLogger()
        self.logger_handler.debug(self.__repr__())

    def run(self):
        tmp_dict = defaultdict(lambda: [])
        self.logger_handler.debug(f"Worker PID: {self.pid} started")
        for bytes_window_start in range(len(self.fasta_bytes) - self.index_len + 1):
            # print(self.start_pos+bytes_window_start)
            bytes_window = self.fasta_bytes[bytes_window_start:bytes_window_start + self.index_len]
            if _is_low_complexity(bytes_window):
                continue

            tmp_dict[bytes_window].append(
                (self.chromosome_name, self.strand, self.start_pos + bytes_window_start)
            )
        dump_target = os.path.join(self.tmp_dir, str(self.pid) + ".pkl")
        logger_handler.debug(f"Worker PID: {self.pid} dumping to {dump_target}...")
        pickle.dump(dict(tmp_dict), open(dump_target, 'wb'))
        logger_handler.debug(f"Worker PID: {self.pid} FIN")
        os.kill(self.pid, 9)  # Suicide

    def __repr__(self):
        try:
            return f"Worker {self.chromosome_name} [{self.strand}] {self.start_pos} -> {self.start_pos + len(self.fasta_bytes)}"
        except AttributeError:
            return "Worker under construction"


# class _MergeWorkerThread(threading.Thread):
#     """
#     Multi-thread merger. Accelerate the process by asynchronously load the dump.
#     """
#     def __init__(self, filename: str, out_dict: dict, merge_mutex: threading.Lock):
#         super().__init__()
#         self.filename = filename
#         self.out_dict = out_dict
#         self.merge_mutex = merge_mutex
#
#     def run(self):
#         logger_handler.debug(f"Merging from {self.filename}..")
#         tmpd = pickle.load(open(self.filename, 'rb'))
#         os.remove(self.filename)
#         with self.merge_mutex:
#             self.out_dict = _blast_index_merge(self.out_dict, tmpd)
#         pass
#


class BlastIndex:
    def __init__(self, pool_len: Optional[int] = multiprocessing.cpu_count()):
        logger_handler.info("BlastIndex initializing...")
        self._index = {}
        """The real index"""

        self._pool_len = pool_len
        """Internal thread number"""

        self._fasta_obj = None
        """Internal Fasta object"""

        self.index_len = None

    def create_index(self, fasta_filename: str, index_len: Optional[int] = DEFAULT_INDEX_LEN):
        """
        Wrapper function to (re)-create the index.

        :param fasta_filename: Filename of input FASTA
        :param index_len: The length (k of the k-mer) of the index.
        """
        logger_handler.info(f"Creating index from {fasta_filename}...")
        self._fasta_obj = in_memory_fasta.Fasta(fasta_filename)

        self.index_len = index_len
        self._index = {}

        pool = parallel_helper.MultiProcessingJobQueue(pool_name="Creating chunk-level Index", pool_size=self._pool_len,
                                                       with_tqdm=True)
        """The pool of process"""

        tmp_dir = tempfile.mkdtemp()
        """Temporary directory for index files"""

        for chromosome_name in self._fasta_obj.chromosomes:
            self._submit_create_index_job_for_chromosome(chromosome_name=chromosome_name, tmp_dir=tmp_dir, pool=pool)

        pool.start()
        pool.join()

        futures=[]
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=self._pool_len)
        logger_handler.info("Submitting unpickle jobs...")
        for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
            futures.append(executor.submit(pickle.load,open(filename,'rb')))
            # tmpdl = list(executor.map(lambda filename: (),)))
        logger_handler.info("Unpickling...")
        executor.shutdown()
        for tmpd in tqdm.tqdm(desc="Merging",iterable=futures):
            self._index=_blast_index_merge(self._index,tmpd.result())

    def show_index(self) -> str:
        rets = ""
        for k, v in self._index.items():
            rets += f"{k}: {v}\n"
        return rets

    def write_index(self, dest_filename: str):
        """Write index to a file."""
        logger_handler.info(f"Writing index to {dest_filename}...")
        compressed_pickle.dump(self._index, dest_filename)
        logger_handler.info(f"Writing index to {dest_filename} FIN")

    def read_index(self, from_filename: str):
        """Read index from a file."""
        logger_handler.info(f"Reading index from {from_filename}...")
        self._index = compressed_pickle.load_with_tqdm(from_filename)
        self.index_len = len(list(self._index.keys())[0])
        logger_handler.info(f"Reading index from {from_filename} FIN")

    def get_stats(self):
        """Print statistics of this index"""
        print(f"Index length: {self.index_len}")
        value_len_dir = defaultdict(lambda: 0)
        for v in self._index.values():
            v_len = len(v)
            value_len_dir[v_len] += 1
        value_len_sum = sum(value_len_dir.values())
        print(f"Value len distribution: {value_len_dir}, sum {value_len_sum}, mean {value_len_sum / len(self._index)}")

    def _submit_create_index_job_for_chromosome(
            self, chromosome_name: str, tmp_dir: str, pool: parallel_helper.MultiProcessingJobQueue
    ) -> None:
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
            input_fasta_bytes = fasta_bytes[chunk_len * i:chunk_len * (i + 1) + self.index_len - 1]
            new_worker = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                             strand=True,
                                             fasta_bytes=input_fasta_bytes,
                                             index_len=self.index_len,
                                             start_pos=chunk_len * i,
                                             tmp_dir=tmp_dir)
            pool.append(new_worker)
            new_worker_rc = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                                strand=False,
                                                fasta_bytes=input_fasta_bytes,
                                                index_len=self.index_len,
                                                start_pos=chunk_len * i,
                                                tmp_dir=tmp_dir)
            pool.append(new_worker_rc)

    def valid_index(self):
        """Validity the index. Should NOT be called by outer process."""
        if self._fasta_obj is None:
            return True
        total_item = sum([len(a) for a in self._index.values()])

        with tqdm.tqdm(desc="Validating index", total=total_item) as pbar:
            for k, v in self._index.items():
                for location in v:
                    assert self._fasta_obj.get(chromosome_name=location[0], start=location[2],
                                               end=location[2] + self.index_len, strand=location[1]) == k
                    pbar.update(1)


if __name__ == '__main__':
    bi = BlastIndex()
    bi.create_index('test/test.fasta')
    # print(bi.show_index())
    bi.valid_index()
    bi.get_stats()
    # bi.write_index('test/test.pkl.xz')
    # bj = BlastIndex()
    # bj.read_index('test/test.pkl.xz')
