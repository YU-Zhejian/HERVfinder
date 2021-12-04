import gc
import itertools
import logging
import multiprocessing
import queue
import time
from collections import defaultdict
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import tqdm

from herv_finder.blast import indexer, blast_anchors_type, \
    blast_merged_anchors_type, blast_merged_anchor_type, \
    blast_extended_anchor_type, blast_extended_anchors_type
from herv_finder.utils import in_memory_fasta

EXTENDER_BATCH_SIZE = 10000
"""After generation of anchors of this number, submit them to the queue."""


def merge_adjacent_anchors(raw_anchors: blast_anchors_type, word_len: int) -> blast_merged_anchors_type:
    """
    Merge adjacent anchors.
    """
    anchor_dict = defaultdict(list)
    """Dict[needle_chromosome, needle_strand, haystack_chromosome, haystack_strand], 
    List[needle_start, haystack_start]]"""
    for anchor in raw_anchors:
        anchor_dict[(anchor[0][0], anchor[0][1], anchor[1][0], anchor[1][1])].append((anchor[0][2], anchor[1][2]))
    for k, v in anchor_dict.items():
        while len(v) > 0:
            tmp_merge_base = v.pop()
            merge_base_f = (tmp_merge_base[0], tmp_merge_base[1])
            merge_base_b = (tmp_merge_base[0], tmp_merge_base[1])
            extend_length = word_len
            # Forward extend
            while True:
                tmp_extend_anchor = (merge_base_f[0] + 1, merge_base_f[1] + 1)
                if tmp_extend_anchor in v:
                    v.remove(tmp_extend_anchor)
                    merge_base_f = (merge_base_f[0] + 1, merge_base_f[1] + 1)
                    extend_length += 1
                else:
                    break
            # Reverse extend
            while True:
                tmp_extend_anchor = (merge_base_b[0] - 1, merge_base_b[1] - 1)
                if tmp_extend_anchor in v:
                    v.remove(tmp_extend_anchor)
                    merge_base_b = (merge_base_b[0] - 1, merge_base_b[1] - 1)
                    extend_length += 1
                else:
                    break
            yield (k[0], k[1], merge_base_b[0]), (k[2], k[3], merge_base_b[1]), extend_length


def blast_score(b1: bytes, b2: bytes) -> int:
    l_b1 = len(b1)
    l_b2 = len(b2)
    scores = np.zeros((l_b1, l_b2), dtype=int)
    for i in range(l_b1):
        for j in range(l_b2):
            if i == 0:
                scores[i][j] = j
            elif j == 0:
                scores[i][j] = i
            elif b1[i - 1] == b2[j - 1]:
                scores[i][j] = 2 + scores[i - 1][j - 1]
            else:
                scores[i][j] = max(
                    scores[i - 1][j - 1],
                    scores[i][j - 1],
                    scores[i - 1][j]
                ) - 1
    return scores[l_b1 - 1][l_b2 - 1]


def visualize_anchors(anchors: blast_merged_anchors_type,
                      needle_chromosome_name: str,
                      haystack_chromosome_name: str,
                      needle_strand: bool = True, haystack_strand: bool = True
                      ):
    x = []
    y = []
    anchors = list(anchors)
    for anchor in tqdm.tqdm(desc="Preparing plot", iterable=anchors):
        if anchor[0][0] == needle_chromosome_name and anchor[1][0] == haystack_chromosome_name and \
                anchor[0][1] == needle_strand and anchor[1][1] == haystack_strand:
            x_start = anchor[0][2]
            y_start = anchor[1][2]
            for i in range(anchor[2]):
                x.append(x_start + i)
                y.append(y_start + i)

    plt.scatter(x, y, s=0.01)
    plt.show()
    plt.hist(x, bins=50)
    plt.show()
    plt.hist(y, bins=50)
    plt.show()


class _ExtendWorkerProcess(multiprocessing.Process):
    def __init__(self, anchors_to_extend: blast_merged_anchors_type,
                 needle_fasta_filename: str,
                 haystack_fasta_name: str,
                 output_queue: multiprocessing.Queue
                 ):
        super().__init__()
        self.anchors_to_extend = anchors_to_extend
        self.needle_fasta = in_memory_fasta.Fasta(needle_fasta_filename)
        self.haystack_fasta = in_memory_fasta.Fasta(haystack_fasta_name)
        self.output_queue = output_queue
        self.logger_handler = logging.getLogger("Worker")
        self.logger_handler.debug(self.__repr__())

    def extend(self, merged_anchor: blast_merged_anchor_type) -> blast_extended_anchor_type:
        # print(merged_anchor)
        needle_len = merged_anchor[2]
        needle_chromosome, needle_strand, needle_start = merged_anchor[0]
        needle_end = needle_start + needle_len
        haystack_chromosome, haystack_strand, haystack_start = merged_anchor[1]
        haystack_end = haystack_start + needle_len
        max_score = needle_len * 2
        max_coordinate = (needle_start, needle_end, haystack_start, haystack_end)
        needle_full_length = self.needle_fasta.get_chromosome_length(needle_chromosome)
        haystack_full_length = self.haystack_fasta.get_chromosome_length(haystack_chromosome)
        needle_starts = range(max(0, needle_start - needle_len), needle_start)
        needle_ends = range(needle_end, min(needle_full_length, needle_end + needle_len))
        haystack_starts = range(max(0, haystack_start - needle_len), haystack_start)
        haystack_ends = range(haystack_end, min(haystack_full_length, haystack_end + needle_len))
        # print(needle_end, min(needle_full_length, needle_end + needle_len))
        # print(list(needle_starts),list(needle_ends),list(haystack_starts),list(haystack_ends))
        for available_coordinate in itertools.product(needle_starts, needle_ends, haystack_starts, haystack_ends):
            # print(available_coordinate)
            needle_bytes = self.needle_fasta.get(needle_chromosome, available_coordinate[0],
                                                 available_coordinate[1], needle_strand)
            haystack_bytes = self.haystack_fasta.get(haystack_chromosome, available_coordinate[2],
                                                     available_coordinate[3], haystack_strand)
            score = blast_score(needle_bytes, haystack_bytes)
            # print(score)
            if score > max_score:
                max_score = score
                max_coordinate = available_coordinate
        return (
            (needle_chromosome, needle_strand, max_coordinate[0], max_coordinate[1] - max_coordinate[0]),
            (haystack_chromosome, haystack_strand, max_coordinate[2], max_coordinate[3] - max_coordinate[2]),
            max_score
        )

    def run(self):
        for alignment in self.anchors_to_extend:
            self.output_queue.put(self.extend(alignment))
        self.logger_handler.info("FIN")

    def __repr__(self):
        try:
            return f"Worker {self.pid}"
        except AttributeError:
            return "Worker under construction"


class BlastIndexSearcher:
    """A memory-efficient BLAST index searcher"""

    def __init__(self,
                 needle_index: indexer.InMemorySimpleBlastIndex,
                 haystack_index: indexer.BlastIndex,
                 pool_len: Optional[int] = multiprocessing.cpu_count(),
                 extend_batch_size:Optional[int] = EXTENDER_BATCH_SIZE
                 ):
        self._pool_len = pool_len
        """Internal thread number"""
        self.haystack_index = haystack_index
        self.needle_index = needle_index
        self.extend_batch_size=extend_batch_size

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
                    yield needle_location, haystack_location

    def extend(self, merged_anchors: blast_merged_anchors_type) -> blast_extended_anchors_type:
        """"""
        this_manager = multiprocessing.Manager()
        output_queue = this_manager.Queue()
        this_batch=[]
        for merged_anchor in merged_anchors:
            if len(this_batch)<self.extend_batch_size:
                this_batch
            else:
                this_batch=[]
        if len(this_batch) > 0:
            this_batch

    def extend_single_thread(self, merged_anchors: blast_merged_anchors_type) -> blast_extended_anchors_type:
        merged_anchors = list(merged_anchors)
        this_manager = multiprocessing.Manager()
        output_queue = this_manager.Queue()
        ewp = _ExtendWorkerProcess(merged_anchors, self.needle_index.fasta_filename,
                                   self.haystack_index.fasta_filename, output_queue)
        ewp.start()
        with tqdm.tqdm(total=len(merged_anchors)) as pbar:
            while ewp.is_alive():
                try:
                    print(ewp.exitcode, ewp.is_alive())
                    yield output_queue.get(timeout=0.1)
                    pbar.update(1)
                    # print(get_item)
                except (TimeoutError, queue.Empty):
                    pass
        ewp.join()
        ewp.close()
        this_manager.shutdown()
        gc.collect()


def _test_tiny():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/tiny_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/tiny")
    haystack_index.attach_fasta("test/tiny.fasta")
    haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index)
    for anchor in searcher.generate_raw_anchors():
        print(anchor[0][2], anchor[1][2])
    for merged_anchor in merge_adjacent_anchors(searcher.generate_raw_anchors(), word_len=11):
        print(merged_anchor)


def _test_on_e_coli():
    needle_index = indexer.InMemorySimpleBlastIndex(word_len=5, prefix_len=1)
    needle_index.attach_fasta("test/e_coli_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/e_coli", word_len=5, prefix_len=1)
    haystack_index.attach_fasta("test/e_coli.fasta")
    # haystack_index.create_index()

    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index)
    raw_anchors = searcher.generate_raw_anchors()
    # print(list(raw_anchors))
    merged_anchors = merge_adjacent_anchors(raw_anchors, word_len=5)
    # searcher.visualize_anchors(merged_anchors,
    #                            'needle', 'NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome')
    for item in searcher.extend_single_thread(itertools.islice(merged_anchors, 100)):
        print(item)


def _test_on_test():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/herv.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/test")
    haystack_index.attach_fasta("test/test.fasta")
    # haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index)
    # for anchor in searcher.generate_raw_anchors():
    #     print(anchor[0][2],anchor[1][2])
    visualize_anchors(merge_adjacent_anchors(searcher.generate_raw_anchors(), 11), 'DF0000558.4 LTR5_Hs', 'chr21')


if __name__ == "__main__":
    _test_on_e_coli()
