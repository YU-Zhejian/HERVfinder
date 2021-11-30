from collections import defaultdict
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

blast_merged_anchor_type=Tuple[blast_index_location_type, blast_index_location_type, int]
blast_merged_anchors_type = Iterable[blast_merged_anchor_type]

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
    return b'N' in fasta_bytes

def merge_adjacent_anchors(raw_anchors:blast_anchors_type) -> blast_merged_anchors_type:
    """
    Merge adjacent anchors.
    """
    anchor_dict = defaultdict(list)
    """Dict[needle_chromosome, needle_strand, haystack_chromosome, haystack_strand], 
    List[needle_start, haystack_start]]"""
    for anchor in raw_anchors:
        anchor_dict[(anchor[0][0],anchor[0][1],anchor[1][0],anchor[1][1])].append((anchor[0][2],anchor[1][2]))
    for k, v in anchor_dict.items():
        while len(v) > 0:
            tmp_merge_base = v.pop()
            merge_base_f = (tmp_merge_base[0], tmp_merge_base[1])
            merge_base_b = (tmp_merge_base[0], tmp_merge_base[1])
            extend_length=11
            # Forward extend
            while True:
                tmp_extend_anchor = (merge_base_f[0] + 1, merge_base_f[1] + 1)
                if tmp_extend_anchor in v:
                    v.remove(tmp_extend_anchor)
                    merge_base_f=(merge_base_f[0] + 1, merge_base_f[1] + 1)
                    extend_length+=1
                else:
                    break
            # Reverse extend
            while True:
                tmp_extend_anchor = (merge_base_b[0] - 1, merge_base_b[1] - 1)
                if tmp_extend_anchor in v:
                    v.remove(tmp_extend_anchor)
                    merge_base_b=(merge_base_b[0] - 1, merge_base_b[1] - 1)
                    extend_length+=1
                else:
                    break
            yield ((k[0], k[1], merge_base_b[0]), (k[2], k[3], merge_base_b[1]), extend_length)
