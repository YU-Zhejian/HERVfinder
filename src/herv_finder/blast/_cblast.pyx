from typing import Tuple

from libc.stdlib cimport malloc, free
from libc.string cimport strlen

blast_index_location_type = Tuple[str, bool, int]

def strand_to_str(strand: bool) -> str:
    """Convert strand information to string"""
    if strand:
        return '+'
    else:
        return '-'

def _editing_distance_blast_score(b1: bytes, b2: bytes) -> int:
    """A complicated blast score by editing distance"""
    return __editing_distance_blast_score(b1, b2)

def _hamming_distance_blast_score(b1: bytes, b2: bytes) -> int:
    """A very simple blast score by hamming distance"""
    return __hamming_distance_blast_score(b1, b2)

cdef int __hamming_distance_blast_score(char * b1, char * b2):
    cdef int l_b1 = strlen(b1)
    cdef int score = 2 * l_b1
    cdef int i
    for i in range(0, l_b1, 1):
        if b1[i] != b2[i]:
            score -= 1
    return score

cdef int __editing_distance_blast_score(char * b1, char * b2):
    cdef int l_b1 = strlen(b1)
    cdef int l_b2 = strlen(b2)
    cdef int * scores
    scores = <int *> malloc(l_b2 * l_b1 * sizeof(int))
    # scores = np.zeros((l_b1, l_b2), dtype=int)
    cdef int i, j
    for i in range(0, l_b1, 1):
        for j in range(0, l_b2, 1):
            if i == 0:
                scores[i*l_b1+j] = j
            elif j == 0:
                scores[i*l_b1+j] = i
            elif b1[i - 1] == b2[j - 1]:
                scores[i*l_b1+j] = 2 + scores[(i-1)*l_b1+(j-1)]
            else:
                scores[i*l_b1+j] = max(
                    scores[(i-1)*l_b1+(j-1)],
                    scores[i*l_b1+(j-1)],
                    scores[(i-1)*l_b1+j]
                ) - 1
    cdef int reti = scores[(i-1)*l_b1+(j-1)]
    free(scores)
    return reti


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
