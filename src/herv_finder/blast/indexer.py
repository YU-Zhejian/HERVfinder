"""A general-purposed multi-processed BLAST index creator"""
import glob
import logging
import math
import multiprocessing
import os
import tempfile
from typing import Tuple, Dict, Optional, Iterator, List

logging.basicConfig(level=logging.DEBUG)
logger_handler = logging.getLogger()

from herv_finder.utils.compressed_pickle import dump, load
from herv_finder.utils.in_memory_fasta import Fasta, get_reversed_complementary

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
        dump_target = os.path.join(self.tmp_dir, str(self.pid) + ".pkl.xz")
        logger_handler.debug(f"Worker PID: {self.pid} dumping to {dump_target}..")
        dump(tmp_dict, dump_target)
        logger_handler.debug(f"Worker PID: {self.pid} FIN")
        os.kill(self.pid, 9)  # Suicide


class BlastIndex():
    def __init__(self, pool_len: Optional[int] = os.cpu_count()):
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
        for chromosome_name in self._fasta_obj.chromosomes:
            self._index.update(self._create_index_for_chromosome(chromosome_name=chromosome_name))

    def write_index(self, dest_filename: str):
        dump(self._index, dest_filename)

    def read_index(self, from_filename: str):

        self._index = load(from_filename)
        self.index_length=len(self._index.keys()[0])

    def get_stats(self):
        print(f"Index length: {self.index_length}")
        value_len_dir={}
        for v in self._index.values():
            v_len = len(v)
            if v_len not in value_len_dir:
                value_len_dir[v_len] =0
            value_len_dir[v_len]+=1
        value_len_sum = sum(value_len_dir.values())
        print(f"value len distribution: {value_len_dir}, sum {value_len_sum}, mean {value_len_sum/len(value_len_dir)}")

    def _create_index_for_chromosome(self, chromosome_name: str) -> blast_index_type:
        tmp_dir = tempfile.mkdtemp()
        """Temporary directory for index files"""

        logger_handler.info(f"Creating index for chromosome {chromosome_name}")
        """Create index for a particular chromosome."""

        # The forward strand
        fasta_bytes = self._fasta_obj.get(chromosome_name)
        """The chromosome in bytes"""

        rc_fasta_bytes = get_reversed_complementary(fasta_bytes)
        """Reverse-complementary of fasta_bytes"""

        total_len = len(fasta_bytes)
        """Total length of this chromosome"""

        chunk_len = math.ceil(total_len / self._pool_len)
        """Length of each chunk"""

        # Split the job
        pool = []
        """The pool of process"""

        # Start the jobs
        for i in range(self._pool_len):
            new_worker = _WorkerProcess(chromosome_name=chromosome_name,
                                        strand=True,
                                        fasta_bytes=fasta_bytes[chunk_len * i:chunk_len * (i + 1) + self.index_length],
                                        index_length=self.index_length,
                                        start_pos=chunk_len * i,
                                        length=chunk_len + self.index_length,
                                        tmp_dir=tmp_dir)
            new_worker.start()
            pool.append(new_worker)
        for process in pool:
            process.join()
        for i in range(self._pool_len):
            new_worker_rc = _WorkerProcess(chromosome_name=chromosome_name,
                                           strand=False,
                                           fasta_bytes=rc_fasta_bytes[
                                                       chunk_len * i:chunk_len * (i + 1) + self.index_length],
                                           index_length=self.index_length,
                                           start_pos=chunk_len * i,
                                           length=chunk_len + self.index_length,
                                           tmp_dir=tmp_dir)
            new_worker_rc.start()
            pool.append(new_worker_rc)
        for process in pool:
            process.join()
        retd = {}
        for filename in glob.glob(os.path.join(tmp_dir, "*.pkl.xz")):
            logger_handler.debug(f"Loading from {filename}..")
            tmpd = load(filename)
            os.remove(filename)
            # print(tmpd.__str__()[1:20])
            retd.update(tmpd)
        os.rmdir(tmp_dir)
        return retd


if __name__ == '__main__':
    bi = BlastIndex()
    bi.create_index('test/sequence.fasta')
    bi.get_stats()
    bi.write_index('test/test.pkl.xz')
