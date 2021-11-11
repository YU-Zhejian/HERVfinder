"""A general-purposed multi-processed BLAST index creator"""
import glob
import logging
import multiprocessing
import os
import pickle
import tempfile
import threading
from typing import Tuple, Dict, Optional, Iterator, List

import tqdm

from herv_finder.utils import multiprocess_helper
from herv_finder.utils import compressed_pickle
from herv_finder.utils.in_memory_fasta import Fasta, get_reversed_complementary

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()

blast_index_type = Dict[bytes, List[Tuple[str, bool, int]]]
"""
The type of blast index. Will be a Dictionary with:

index_item -> [chromosome_name, strand, integer]

strand is True for default and False for reversed.
"""


class _WorkerProcess(multiprocessing.Process):
    def __init__(self,
                 chromosome_name: str,
                 strand: bool,
                 fasta_bytes: bytes,
                 index_length: int,
                 start_pos: int,
                 length: int,
                 tmp_dir: str
                 ):
        """
        The worker of index creator. Will create index for both forward and reverse strand.

        :param chromosome_name: Name of chromosome.
        :param fasta_bytes: Bytes needed to be indexed.
        """
        super().__init__()
        logger_handler.debug(f"Worker add: {chromosome_name} [{strand}] {start_pos} -> {start_pos + length}")
        self.length = length
        self.start_pos = start_pos
        self.index_length = index_length
        self.fasta_bytes = fasta_bytes
        self.strand = strand
        self.chromosome_name = chromosome_name
        self.tmp_dir = tmp_dir

    @staticmethod
    def _is_low_complexity(fasta_bytes: bytes) -> bool:
        return b'N' in fasta_bytes or b'a' in fasta_bytes

    def run(self):
        tmp_dict = {}
        logger_handler.debug(f"Worker PID: {self.pid} started")
        for bytes_window_start in range(0, self.length + self.index_length):
            bytes_window = self.fasta_bytes[bytes_window_start:bytes_window_start + self.index_length]
            if self._is_low_complexity(bytes_window):
                continue
            if bytes_window not in tmp_dict:
                tmp_dict[bytes_window] = []
            tmp_dict[bytes_window].append(((self.chromosome_name, self.strand, self.start_pos + bytes_window_start)))
            # print(self.pid, (self.chromosome_name, self.strand, self.start_pos+ bytes_window_start))
        dump_target = os.path.join(self.tmp_dir, str(self.pid) + ".pkl")
        logger_handler.debug(f"Worker PID: {self.pid} dumping to {dump_target}..")
        pickle.dump(tmp_dict, open(dump_target,'wb'))
        logger_handler.debug(f"Worker PID: {self.pid} FIN")
        os.kill(self.pid, 9)  # Suicide

class _MergeWorkerThread(threading.Thread):
    def __init__(self, filename:str, out_dict:dict, merge_mutex:threading.Lock):
        super().__init__()
        self.filename = filename
        self.out_dict = out_dict
        self.merge_mutex=merge_mutex

    def run(self):
        # logger_handler.debug(f"Loading from {self.filename}..")
        tmpd = pickle.load(open(self.filename, 'rb'))
        os.remove(self.filename)
        with self.merge_mutex:
            self.out_dict.update(tmpd)

class BlastIndex():
    def __init__(self, pool_len: Optional[int] = multiprocessing.cpu_count()):
        logger_handler.info("BlastIndex initializing...")
        self._index = {}
        """The real index"""

        self._pool_len = pool_len
        """Internal thread number"""

    @staticmethod
    def _merge(unmerged_index: Iterator[blast_index_type]) -> blast_index_type:
        """ To merge all unmerged index (passed in as an Iterator) to one. """
        pass

    def create_index(self, fasta_filename: str, index_length: Optional[int] = 11):
        """
        Wrapper function to (re)-create the index.

        :param fasta_filename: Filename of input FASTA
        :param index_length: Length of the index.
        """
        logger_handler.info(f"Creating index from {fasta_filename}...")
        self._fasta_obj = Fasta(fasta_filename)
        """Internal Fasta object"""

        self.index_length = index_length
        self._index = {}

        pool = multiprocess_helper.MPJobQueue(pool_size=self._pool_len, with_tqdm=True)
        """The pool of process"""

        tmp_dir = tempfile.mkdtemp()
        """Temporary directory for index files"""

        for chromosome_name in self._fasta_obj.chromosomes:
            self._submit_create_index_job_for_chromosome(chromosome_name=chromosome_name, tmp_dir=tmp_dir, pool=pool)

        pool.start()
        pool.join()

        logger_handler.info("Merging indicies...")
        merge_pool=[]
        merge_mutex = threading.Lock()

        for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
            merge_pool.append(_MergeWorkerThread(filename=filename, out_dict=self._index, merge_mutex = merge_mutex))
            merge_pool[-1].start()
        for thread in tqdm.tqdm(iterable=merge_pool):
            thread.join()
        os.rmdir(tmp_dir)


    def write_index(self, dest_filename: str):
        logger_handler.info(f"Writing index to {dest_filename}...")
        compressed_pickle.dump(self._index, dest_filename)
        logger_handler.info(f"Writing index to {dest_filename} FIN")

    def read_index(self, from_filename: str):
        logger_handler.info(f"Reading index from {from_filename}...")
        self._index = compressed_pickle.load_with_tqdm(from_filename)
        self.index_length = len(list(self._index.keys())[0])
        logger_handler.info(f"Reading index from {from_filename} FIN")

    def get_stats(self):
        print(f"Index length: {self.index_length}")
        value_len_dir = {}
        for v in self._index.values():
            v_len = len(v)
            if v_len not in value_len_dir:
                value_len_dir[v_len] = 0
            value_len_dir[v_len] += 1
        value_len_sum = sum(value_len_dir.values())
        print(f"value len distribution: {value_len_dir}, sum {value_len_sum}, mean {value_len_sum / len(self._index)}")

    def _submit_create_index_job_for_chromosome(self, chromosome_name: str,tmp_dir:str, pool:multiprocess_helper.MPJobQueue) -> None:
        """Submit index creation job for a particular chromosome."""

        # The forward strand
        fasta_bytes = self._fasta_obj.get(chromosome_name)
        """The chromosome in bytes"""

        rc_fasta_bytes = get_reversed_complementary(fasta_bytes)
        """Reverse-complementary of fasta_bytes"""

        total_len = len(fasta_bytes)
        """Total length of this chromosome"""

        chunk_len = 100000
        """Length of each chunk"""

        # Submit the jobs
        for i in range(total_len // chunk_len + 1):
            new_worker = _WorkerProcess(chromosome_name=chromosome_name,
                                        strand=True,
                                        fasta_bytes=fasta_bytes[chunk_len * i:chunk_len * (i + 1) + self.index_length],
                                        index_length=self.index_length,
                                        start_pos=chunk_len * i,
                                        length=chunk_len + self.index_length,
                                        tmp_dir=tmp_dir)
            pool.append(new_worker)
        for i in range(self._pool_len):
            new_worker_rc = _WorkerProcess(chromosome_name=chromosome_name,
                                           strand=False,
                                           fasta_bytes=rc_fasta_bytes[
                                                       chunk_len * i:chunk_len * (i + 1) + self.index_length],
                                           index_length=self.index_length,
                                           start_pos=chunk_len * i,
                                           length=chunk_len + self.index_length,
                                           tmp_dir=tmp_dir)
            pool.append(new_worker_rc)



if __name__ == '__main__':
    bi = BlastIndex()
    bi.create_index('test/sequence.fasta')
    bi.get_stats()
    # bi.write_index('test/test.pkl.xz')
    # bj = BlastIndex()
    # bj.read_index('test/test.pkl.xz')
