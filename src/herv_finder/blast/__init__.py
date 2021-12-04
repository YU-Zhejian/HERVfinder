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
"""
[[chromosome_name, strand, offset],[chromosome_name, strand, offset]]
"""

blast_anchors_type = Iterable[blast_anchor_type]
"""
[[[chromosome_name, strand, offset],[chromosome_name, strand, offset]]]
"""

blast_extended_anchor_type = Tuple[Tuple[str, bool, int, int], Tuple[str, bool, int, int], int]
"""
[[chromosome_name, strand, offset,len],[chromosome_name, strand, offset,len],score]
"""

blast_extended_anchors_type = Iterable[blast_extended_anchor_type]
"""
[[[chromosome_name, strand, offset,len],[chromosome_name, strand, offset,len],score]]
"""

DEFAULT_WORD_LEN = 11
DEFAULT_CHUNK_LEN = 2000000
DEFAULT_PREFIX_LEN = 3
