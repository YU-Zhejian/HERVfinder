import multiprocessing
from typing import Iterable, Tuple, Optional

from herv_finder.blast import blast_index_location_type, blast_simple_index_type, indexer

blast_anchor_type = Tuple[blast_index_location_type, blast_index_location_type]

blast_anchors_type = Iterable[blast_anchor_type]


class _ExtendWorkerProcess(multiprocessing.Process):
    pass


class BlastIndexSearcher:
    """A memory-efficient BLAST index searcher"""

    def __init__(self,
                 needle_index: indexer.InMemorySimpleBlastIndex,
                 basename: str,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 ):
        self.index_basename = basename
        self._pool_len = pool_len
        """Internal thread number"""
        self._haystack_index = indexer.BlastIndex(
            basename=basename,
            pool_len=pool_len
        )
        self.needle_index = needle_index

    def _load_segment(self, prefix_fasta_bytes: bytes) -> blast_simple_index_type:
        """"""
        pass

    def search(self, input_positions: blast_simple_index_type) -> blast_anchors_type:
        """Search a blast_simple_index_type against another index"""
        pass

    def extend(self, search_input: blast_anchors_type) -> blast_anchors_type:
        """"""
        pass
