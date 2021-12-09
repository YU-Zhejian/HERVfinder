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

DEFAULT_WORD_LEN = 11
DEFAULT_CHUNK_LEN = 2000000
DEFAULT_PREFIX_LEN = 3
DEFAULT_SCORE_CUTOFF_SLOPE = 1.5
DEFAULT_SCORE_CUTOFF_INTERSECT = 0.0
DEFAULT_EXTENDER_BATCH_SIZE = 100


def strand_to_str(strand: bool) -> str:
    """Convert strand information to string"""
    if strand:
        return '+'
    else:
        return '-'


class BlastAnchor:
    __slots__ = [
        "haystack_start",
        "haystack_end",
        "haystack_chromosome",
        "haystack_strand",
        "needle_start",
        "needle_end",
        "needle_chromosome",
        "score",
        "strand"
    ]
    haystack_start: int
    haystack_end: int
    haystack_chromosome: str
    haystack_strand: bool
    needle_start: int
    needle_end: int
    needle_chromosome: str
    score: int
    strand: bool

    def __init__(self):
        self.score = 0

    @classmethod
    def from_params(cls,
                    haystack_start: int = 0,
                    haystack_end: int = 0,
                    haystack_chromosome: str = "chr?",
                    needle_start: int = 0,
                    needle_end: int = 0,
                    needle_chromosome: str = "chr?",
                    score: int = 0,
                    strand: bool = True,
                    ):
        instance = cls()
        instance.haystack_start = haystack_start
        instance.haystack_end = haystack_end
        instance.haystack_chromosome = haystack_chromosome
        instance.needle_start = needle_start
        instance.needle_end = needle_end
        instance.needle_chromosome = needle_chromosome
        instance.score = score
        instance.strand = strand
        return instance

    @classmethod
    def from_location(cls,
                      needle_location: blast_index_location_type,
                      haystack_location: blast_index_location_type,
                      word_len: int):
        instance = cls()
        instance.needle_chromosome = needle_location[0]
        instance.strand = needle_location[1] and haystack_location[1]
        instance.needle_start = needle_location[2]
        instance.needle_end = instance.needle_start + word_len
        instance.haystack_chromosome = haystack_location[0]
        instance.haystack_start = haystack_location[2]
        instance.haystack_end = instance.haystack_start + word_len
        return instance

    def to_gtf(self) -> str:
        return "\t".join((
            self.haystack_chromosome,  # Chromosome
            "HERVfinder",  # Source
            "exon",  # Type
            str(self.haystack_start),  # Start
            str(self.haystack_end),  # End
            str(self.score),  # Score
            strand_to_str(self.strand),  # Strand
            ".",  # Frame
            f"gene_id \"{self.needle_chromosome}\";"  # Suppl
        ))

    @property
    def needle_len(self) -> int:
        return self.needle_end - self.needle_start

    @property
    def haystack_len(self) -> int:
        return self.haystack_end - self.haystack_start

    def __str__(self) -> str:
        return ",".join((
            self.needle_chromosome,
            str(self.needle_start),
            str(self.needle_end),
            self.haystack_chromosome,
            str(self.haystack_start),
            str(self.haystack_end),
            strand_to_str(self.strand),
            str(self.score)
        ))
