from typing import List, Dict, Tuple

blast_index_location_type = Tuple[str, bool, int]
"""
The basic index segment.

[chromosome_name, strand, offset]

strand is True for default and False for reversed.
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
Indexed using leading some bases to save memory.
"""

DEFAULT_WORD_LEN = 11
"""Index length of nucleotides"""

DEFAULT_CHUNK_LEN = 2000000
"""Default chunk length"""

DEFAULT_LEADING_N_BASES = 1
""""""

def blast_index_merge(d1: blast_simple_index_type, d2: blast_simple_index_type):
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



