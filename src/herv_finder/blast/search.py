import itertools
import multiprocessing
from typing import Iterable, Tuple, Optional

from herv_finder.blast import blast_index_location_type, blast_simple_index_type, indexer, blast_anchors_type


class _ExtendWorkerProcess(multiprocessing.Process):
    pass


class BlastIndexSearcher:
    """A memory-efficient BLAST index searcher"""

    def __init__(self,
                 needle_index: indexer.InMemorySimpleBlastIndex,
                 haystack_index: indexer.BlastIndex,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 ):
        self._pool_len = pool_len
        """Internal thread number"""
        self.haystack_index = haystack_index
        self.needle_index = needle_index
        self.sorted_anchors = []

    def generate_raw_anchors(self) -> blast_anchors_type:
        """Generate a series of anchors"""
        haystack_available_prefixes = self.haystack_index.available_prefixes
        needle_words = self.needle_index.loaded_words
        for needle_prefix in self.needle_index.loaded_prefixes:
            if needle_prefix not in haystack_available_prefixes:
                continue
            self.haystack_index.read_index([needle_prefix])
            for needle_word in needle_words:
                if not needle_word.startswith(needle_prefix):
                    continue
                for needle_location, haystack_location in itertools.product(
                        self.needle_index.get(needle_word),
                    self.haystack_index.get(needle_word)
                ):
                    yield (needle_location, haystack_location)



def _test_tiny():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/tiny_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/tiny")
    haystack_index.attach_fasta("test/tiny.fasta")
    haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index,haystack_index=haystack_index)
    print(list(searcher.anchor()))



if __name__ == "__main__":
    _test_tiny()
