"""A general-purposed multi-processed BLAST index creator"""
import glob
import logging
import multiprocessing
import os
import pickle
import tempfile
import uuid
from collections import defaultdict
from typing import Optional, List

import tqdm

from herv_finder.blast import DEFAULT_WORD_LEN, DEFAULT_CHUNK_LEN, \
    blast_index_merge, is_low_complexity, DEFAULT_LEADING_N_BASES
from herv_finder.utils import compressed_pickle
from herv_finder.utils import in_memory_fasta
from herv_finder.utils import parallel_helper


class _IndexWorkerProcess(multiprocessing.Process):
    def __init__(self,
                 chromosome_name: str,
                 fasta_bytes: bytes,
                 word_len: int,
                 start_pos: int,
                 tmp_dir: str,

                 ):
        """
        The worker of index creator. Will create index for both forward and reverse strand.

        :param chromosome_name: Name of chromosome.
        :param fasta_bytes: Bytes needed to be indexed.
        :param word_len: The length (k of the k-mer) of the index.
        :param start_pos: 0-based offset of the fasta_bytes.
        """
        super().__init__()
        self.start_pos = start_pos
        self.word_len = word_len
        self.fasta_bytes = fasta_bytes
        self.chromosome_name = chromosome_name
        self.tmp_dir = tmp_dir
        self.logger_handler = logging.getLogger("Worker")
        self.logger_handler.debug(self.__repr__())
        self.tmp_dict = None

    def run(self):
        self.tmp_dict = defaultdict(lambda: [])
        self.logger_handler.debug(f"Worker PID: {self.pid} started")
        for bytes_window_start in range(len(self.fasta_bytes) - self.word_len + 1):
            # print(self.start_pos+bytes_window_start)
            bytes_window = self.fasta_bytes[bytes_window_start:bytes_window_start + self.word_len]
            if not is_low_complexity(bytes_window):
                self.tmp_dict[bytes_window].append(
                    (self.chromosome_name, True, self.start_pos + bytes_window_start)
                )
            # bytes_window_rc = in_memory_fasta.get_reversed_complementary(bytes_window)
            # if not is_low_complexity(bytes_window_rc):
            #     self.tmp_dict[bytes_window_rc].append(
            #         (self.chromosome_name, False, self.start_pos + bytes_window_start)
            #     )
        self.logger_handler.debug(f"Worker PID: {self.pid} start transmitting data")
        total_items_get = uuid.uuid4()
        for two_mer in in_memory_fasta.get_all_kmers(DEFAULT_LEADING_N_BASES, bases=b'AGCTN'):
            two_mer_str = str(two_mer, encoding='UTF-8')
            index_started_by_two_mer = {}
            for k, v in self.tmp_dict.items():
                if k.startswith(two_mer):
                    index_started_by_two_mer[k] = v
            if len(index_started_by_two_mer) > 0:
                output_filename = f"{two_mer_str}-{total_items_get}.pkl"
                pickle.dump(index_started_by_two_mer,
                            open(os.path.join(self.tmp_dir, output_filename), "wb"))
                del index_started_by_two_mer
        del self.tmp_dict
        self.logger_handler.debug(f"Worker PID: {self.pid} FIN")

    def __repr__(self):
        try:
            return f"Worker {self.chromosome_name} " + \
                   f"{self.start_pos} -> {self.start_pos + len(self.fasta_bytes)}"
        except AttributeError:
            return "Worker under construction"


class _BaseBlastIndex:
    def __init__(self,
                 word_len: Optional[int] = DEFAULT_WORD_LEN,
                 ):
        self.word_len = word_len

        self._indices = {}
        """
        A dict of re-split indices. Format:

        Dict[leading-bytes: blast_simple_index_type]
        """
        self.logger_handler = logging.getLogger()
        self._fasta_obj = None
        """Internal FASTA object"""

        self.total_index_count = 0
        """Number of entries inside the indices. Can be duplicated."""

    def __str__(self) -> str:
        rets = ""
        for two_mer, index in self._indices.items():
            rets += str(two_mer, encoding='UTF-8') + '\n'
            for k, v in index.items():
                rets += f"{k}: {v}\n"
            rets += '\n'
        return rets

    def __repr__(self):
        try:
            return f"Index of {self._fasta_obj.filename}"
        except AttributeError:
            return "Index of unknown FASTA"

    def detach_fasta(self):
        """Detach and free the atached FASTA object"""
        if self._fasta_obj is not None:
            del self._fasta_obj
            self._fasta_obj = None

    def attach_fasta(self, fasta_filename: str):
        """Attach a FASTA object."""
        self._fasta_obj = in_memory_fasta.Fasta(fasta_filename)

    def get_stats(self) -> str:
        """Print statistics of this index"""
        rets=f"Index length: {self.word_len}"
        self._update_total_index_count()
        if self.total_index_count == 0:
            return "None"
        value_len_dir = defaultdict(lambda: 0)
        value_len_sum = 0
        for index in self._indices.values():
            for v in index.values():
                v_len = len(v)
                value_len_dir[v_len] += 1
                value_len_sum += v_len
        rets+= f"Value len distribution: sum {value_len_sum}, mean {value_len_sum / self.total_index_count}"
        rets+= str(dict(value_len_dir))
        return rets

    def _update_total_index_count(self):
        self.total_index_count = 0
        for index in tqdm.tqdm(desc="Updating total index count",
                               iterable=self._indices.values()):
            for v in index.values():
                self.total_index_count += len(v)

    def validate_index(self):
        """Validate the index. Should NOT be called by outer process."""
        if self._fasta_obj is None:
            return

        with tqdm.tqdm(desc="Validating index", total=self.total_index_count) as pbar:
            for index in self._indices.values():
                for k, v in index.items():
                    for location in v:
                        assert self._fasta_obj.get(
                            chromosome_name=location[0],
                            start=location[2],
                            end=location[2] + self.word_len,
                            strand=location[1]) == k
                        pbar.update(1)

    @property
    def keys_full_length(self):
        """All recorded kmers"""
        return [v.keys() for v in self._indices.values()]


class InMemorySimpleBlastIndex(_BaseBlastIndex):
    def __init__(self,
                 word_len: Optional[int] = DEFAULT_WORD_LEN,
                 ):

        super().__init__(word_len=word_len)
        self.logger_handler.info("InMemorySimpleBlastIndex initializing...")

    def create_index(self):
        tmp_dict = defaultdict(lambda: [])
        self._indices = {}
        with tqdm.tqdm(total=self._fasta_obj.total_length) as pbar:
            for chromosome_name in self._fasta_obj.chromosomes:
                """Submit index creation job for a particular chromosome."""
                # The forward strand
                fasta_bytes = self._fasta_obj.get(chromosome_name)
                for bytes_window_start in range(len(fasta_bytes) - self.word_len + 1):
                    # print(self.start_pos+bytes_window_start)
                    bytes_window = fasta_bytes[bytes_window_start:bytes_window_start + self.word_len]
                    if not is_low_complexity(bytes_window):
                        tmp_dict[bytes_window].append(
                            (chromosome_name, True, bytes_window_start)
                        )
                    pbar.update(1)
                    # bytes_window_rc = in_memory_fasta.get_reversed_complementary(bytes_window)
                    # if not is_low_complexity(bytes_window_rc):
                    #     tmp_dict[bytes_window_rc].append(
                    #         (chromosome_name, False, bytes_window_start)
                    #     )
                    # pbar.update(1)
        self.logger_handler.info(f"Re-merging in prefix-level...")
        for two_mer in tqdm.tqdm(desc="Merging...",
                                 iterable=list(in_memory_fasta.get_all_kmers(DEFAULT_LEADING_N_BASES, bases=b'AGCTN'))):
            index_started_by_two_mer = {}
            for k, v in tmp_dict.items():
                if k.startswith(two_mer):
                    index_started_by_two_mer[k] = v
            if len(index_started_by_two_mer) > 0:
                self._indices[two_mer] = index_started_by_two_mer
        self._update_total_index_count()


class BlastIndex(_BaseBlastIndex):
    """The blast index creator and utilities"""

    def __init__(self,
                 basename: str,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 word_len: Optional[int] = DEFAULT_WORD_LEN,
                 chunk_len: Optional[int] = DEFAULT_CHUNK_LEN,
                 ):
        super().__init__(word_len=word_len)

        self.pool_len = pool_len
        """Internal thread number"""

        self.chunk_len = chunk_len

        self.basename = basename
        """The basename of the index"""

    def create_index(self):
        """
        Wrapper function to (re)-create the index.
        """

        self._indices = {}

        pool = parallel_helper.ParallelJobQueue(
            pool_name="Creating chunk-level Index",
            pool_size=self.pool_len,
            with_tqdm=True,
        )
        """The pool of process"""

        tmp_dir = tempfile.mkdtemp()
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
                                    start_pos + self.chunk_len + self.word_len - 1
                                    ]
                new_worker = _IndexWorkerProcess(chromosome_name=chromosome_name,
                                                 fasta_bytes=input_fasta_bytes,
                                                 word_len=self.word_len,
                                                 start_pos=start_pos,
                                                 tmp_dir=tmp_dir)
                pool.append(new_worker)

        pool.start()
        pool.join()
        self.logger_handler.info(f"Chunk-level indexing had finished, Re-merging in prefix-level...")
        with tqdm.tqdm(desc="Merging...",
                       total=len(list(glob.glob(os.path.join(tmp_dir, f"*-*.pkl"))))) as pbar:
            for two_mer in in_memory_fasta.get_all_kmers(DEFAULT_LEADING_N_BASES, bases=b'AGCTN'):
                two_mer_str = str(two_mer, encoding='UTF-8')
                index_started_by_two_mer = {}
                for filename in glob.glob(os.path.join(tmp_dir, f"{two_mer_str}-*.pkl")):
                    self.logger_handler.debug(f"Joining chunk-level index from {filename}")
                    get_item = pickle.load(open(filename, "rb"))
                    blast_index_merge(index_started_by_two_mer, get_item)
                    os.remove(filename)
                    pbar.update(1)
                    del get_item
                if len(index_started_by_two_mer) > 0:
                    output_filename = f"{self.basename}_{two_mer_str}.pkl.xz"
                    self.logger_handler.debug(f"Writing index to {output_filename}...")
                    compressed_pickle.dump(index_started_by_two_mer, output_filename)

        os.rmdir(tmp_dir)
        self.logger_handler.info(f"Index FIN")

    def read_index(self, two_mers: List[bytes] = None):
        """Read index from a file."""
        if two_mers is None:
            two_mers = in_memory_fasta.get_all_kmers(DEFAULT_LEADING_N_BASES, bases=b'AGCTN')
        self.logger_handler.info(f"Reading index from {self.basename}...")
        self._indices = {}
        for two_mer in two_mers:
            kmer_str = str(two_mer, encoding='UTF-8')
            load_filename = f"{self.basename}_{kmer_str}.pkl.xz"
            # print(f"{self.basename}_{kmer_str}.pkl.xz")
            if os.path.exists(load_filename):
                self.logger_handler.info(f"Reading index from {load_filename}...")
                self._indices[two_mer] = compressed_pickle.load_with_tqdm(load_filename)
        try:
            self.word_len = len(list(list(self._indices.values())[0].keys())[0])
        except (IndexError, AttributeError):
            self.word_len = None
        self._update_total_index_count()
        self.logger_handler.info(f"Reading index from {self.basename} FIN")

    def drop_index(self, two_mers: List[bytes] = in_memory_fasta.get_all_kmers(DEFAULT_LEADING_N_BASES, bases=b'AGCTN')):
        for two_mer in two_mers:
            try:
                self._indices.pop(two_mer)
            except KeyError:
                pass
        self._update_total_index_count()


def _base_test(bi: BlastIndex):
    bi.create_index()
    bi.read_index()
    bi.validate_index()
    bi.detach_fasta()
    print(bi.get_stats())


def _test_on_tiny():
    bi = BlastIndex(basename='test/tiny', word_len=4, chunk_len=20)
    bi.attach_fasta('test/tiny.fasta')
    _base_test(bi)
    # print(bi.show_index())
    sbi = InMemorySimpleBlastIndex(word_len=4)
    sbi.attach_fasta('test/tiny.fasta')
    sbi.create_index()
    sbi.validate_index()
    sbi.detach_fasta()
    print(sbi.get_stats())


def _test_on_e_coli():
    bi = BlastIndex(basename='test/e_coli', word_len=4)
    bi.attach_fasta('test/e_coli.fasta')
    _base_test(bi)


def _test_on_test():
    bi = BlastIndex(basename='test/test', pool_len=multiprocessing.cpu_count() // 1)
    bi.attach_fasta('test/test.fasta')
    _base_test(bi)


if __name__ == '__main__':
    for fn in glob.glob('test/*.pkl.xz'):
        os.remove(fn)
    _test_on_tiny()
    # _test_on_e_coli()
    # _test_on_test()
