import multiprocessing
from typing import Iterable, Tuple, Optional

from herv_finder.blast import blast_index_segment_type, blast_index_type

blast_anchor_type = Tuple[blast_index_segment_type, blast_index_segment_type]

blast_anchors_type = Iterable[blast_anchor_type]


class _ExtendWorkerProcess(multiprocessing.Process):
    pass


class BlastIndexSearcher:
    """A memory-efficient BLAST index searcher"""

    def __init__(self,
                 index_basename: str,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 ):
        self.index_basename = index_basename
        self._pool_len = pool_len
        """Internal thread number"""

    def _load_segment(self, prefix_fasta_bytes: bytes) -> blast_index_type:
        pass

    def search(self, input_positions: blast_index_type) -> blast_anchors_type:
        """Search a blast_index_type against another index"""

    def extend(self, search_input: blast_anchors_type) -> blast_anchors_type:
        pass
