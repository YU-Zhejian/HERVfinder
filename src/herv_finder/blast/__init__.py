from typing import List, Dict, Tuple

blast_index_segment_type = Tuple[str, bool, int]
"""
The basic index segment.

[chromosome_name, strand, offset]

strand is True for default and False for reversed.
"""

blast_index_data_type = List[blast_index_segment_type]
"""
The values of `blast_index_type`. Should be a list of `blast_index_segment_type`
"""

blast_index_type = Dict[bytes, blast_index_data_type]
"""
The type of blast index. Will be a Dictionary with:

index_item -> blast_index_data_type
"""

DEFAULT_INDEX_LEN = 11
"""Index length of nucleotides"""

DEFAULT_CHUNK_LEN = 2000000

"""Default chunk length"""

