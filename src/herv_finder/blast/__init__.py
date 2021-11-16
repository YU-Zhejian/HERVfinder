from typing import List, Dict, Tuple, Iterable

blast_index_location_type = Tuple[str, bool, int]
"""
Type of a location. Will cost 64 bytes.

[chromosome_name, strand, offset]

strand is True for default and False for reverse complementary strand.
"""

blast_index_locations_type = List[blast_index_location_type]
"""
The values of `blast_simple_index_type`. Should be a list of `blast_index_location_type`
"""

blast_simple_index_type = Dict[bytes, blast_index_locations_type]
"""
The type of blast index. Will be a Dictionary with:

index_item -> blast_index_locations_type
"""

blast_indices_type = Dict[bytes, blast_simple_index_type]
"""
Indexed using prefix to save memory.
"""

blast_anchor_type = Tuple[blast_index_location_type, blast_index_location_type]

blast_anchors_type = Iterable[blast_anchor_type]


DEFAULT_WORD_LEN = 11
DEFAULT_CHUNK_LEN = 2000000
DEFAULT_PREFIX_LEN = 3


def merge_blast_simple_index(d1: blast_simple_index_type, d2: blast_simple_index_type):
    """Merge d2 to d1"""
    for k, v in d2.items():
        if k in d1:
            d1[k].extend(v)
        else:
            d1[k] = v


def is_low_complexity(fasta_bytes: bytes) -> bool:
    """
    Whether the sequence is of low complexity.

    TODO
    """
    # return b'N' in fasta_bytes or b'a' in fasta_bytes
    _ = fasta_bytes
    return False

def sort_raw_anchors(raw_anchors:blast_anchors_type) -> blast_anchors_type:
    # TODO
    raw_anchor_list=list(raw_anchors)
    pass

def merge_adjacent_anchors(sorted_anchors:blast_anchors_type) -> blast_anchors_type:
    # TODO
    pass
