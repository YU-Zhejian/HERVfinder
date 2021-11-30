import itertools
import multiprocessing
from typing import Iterable, Tuple, Optional

import matplotlib.pyplot as plt
import tqdm

from herv_finder.blast import blast_index_location_type, blast_simple_index_type, indexer, blast_anchors_type, \
    merge_adjacent_anchors, blast_anchor_type, blast_merged_anchors_type
from herv_finder.utils import in_memory_fasta

EXTENDER_BATCH_SIZE = 10000
"""After generation of anchors of this number, submit them to the queue."""

class _ExtendWorkerProcess(multiprocessing.Process):
    def __init__(self,anchors_to_extend:blast_merged_anchors_type,
                 needle_fasta_filename:str,
                 haystack_fasta_name:str,
                 alignment_output:multiprocessing.Queue
                 ):
        super().__init__()
        self.anchors_to_extend = anchors_to_extend
        self.needle_fasta = in_memory_fasta.Fasta(needle_fasta_filename)
        self.haystack_fasta=in_memory_fasta.Fasta(haystack_fasta_name)

    def extend(self, alignment:blast_anchor_type):
        pass

    def run(self):
        retl = [self.extend(alignment) for alignment in self.anchors_to_extend]


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
        for needle_prefix in tqdm.tqdm(iterable=self.needle_index.loaded_prefixes):
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

    def visualize_anchors(self ,
                          anchors:blast_merged_anchors_type,
                          needle_chromosome_name: str,
                          haystack_chromosome_name: str,
                          needle_strand:bool=True, haystack_strand:bool=True
                          ):
        x=[]
        y=[]
        for anchor in tqdm.tqdm(desc="Preparing plot",iterable=anchors):
            if anchor[0][0] == needle_chromosome_name and anchor[1][0] == haystack_chromosome_name and \
                anchor[0][1] == needle_strand and anchor[1][1] ==haystack_strand:
                x_start = anchor[0][2]
                y_start = anchor[1][2]
                for i in range(anchor[2]):
                    x.append(x_start+i)
                    y.append(y_start + i)

        plt.scatter(x,y,s=0.01)
        plt.show()
        plt.hist(x,bins=50)
        plt.show()
        plt.hist(y,bins=50)
        plt.show()


def _test_tiny():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/tiny_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/tiny")
    haystack_index.attach_fasta("test/tiny.fasta")
    haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index,haystack_index=haystack_index)
    for anchor in searcher.generate_raw_anchors():
        print(anchor[0][2],anchor[1][2])
    for merged_anchor in merge_adjacent_anchors(searcher.generate_raw_anchors()):
        print(merged_anchor)

def _test_on_test():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/herv.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/test")
    haystack_index.attach_fasta("test/test.fasta")
    # haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index,haystack_index=haystack_index)
    # for anchor in searcher.generate_raw_anchors():
    #     print(anchor[0][2],anchor[1][2])
    searcher.visualize_anchors(merge_adjacent_anchors(searcher.generate_raw_anchors()),'DF0000558.4 LTR5_Hs','chr21')



if __name__ == "__main__":
    _test_on_test()
