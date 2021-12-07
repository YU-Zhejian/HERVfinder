import gc
import itertools
import logging
import multiprocessing
import os
import pickle
import queue
import tempfile
from collections import defaultdict
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import tqdm

from herv_finder.blast import indexer, blast_anchors_type, \
    blast_extended_anchor_type, blast_extended_anchors_type, blast_anchor_type
from herv_finder.utils import in_memory_fasta, parallel_helper

EXTENDER_BATCH_SIZE = 100
"""After generation of anchors of this number, submit them to the queue."""


def strand_to_str(strand:bool) -> str:
    if strand:
        return '+'
    else:
        return '-'

def extended_anchor_to_gtf(extended_anchor: blast_extended_anchor_type) -> str:
    return f"{extended_anchor[1][0]}\tHERVfinder\texon\t{extended_anchor[1][2]}\t{extended_anchor[1][2]+extended_anchor[1][3]}\t{extended_anchor[2]}\t{strand_to_str(extended_anchor[1][1])}\t.\tgene_id \"{extended_anchor[0][0]}\";"

def simple_blast_scores(b1: bytes, b2: bytes) -> int:
    assert len(b1) == len(b2)
    score = 2 * len(b1)
    for i_b, j_b in zip(b1, b2):
        if i_b != j_b:
            score -= 1
    return score


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


def visualize_anchors(anchors: blast_anchors_type,
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


def merge_adjacent_anchors(
        extended_anchors: blast_extended_anchors_type
) -> blast_extended_anchors_type:
    """
    Merge adjacent anchors.
    """
    anchor_dict = defaultdict(set)
    """Dict[needle_chromosome, needle_strand, haystack_chromosome, haystack_strand], 
    List[needle_start, needle_len, haystack_start, haystack_len]]"""
    for anchor in extended_anchors:
        anchor_dict[(anchor[0][0], anchor[0][1], anchor[1][0], anchor[1][1])].add(
            (anchor[0][2], anchor[0][3], anchor[1][2], anchor[1][3])
        )
    for k, v in anchor_dict.items():
        while len(v) > 0:
            tmp_merge_base = v.pop()
            yield (k[0], k[1], tmp_merge_base[0],tmp_merge_base[1]), (k[2], k[3], tmp_merge_base[2],tmp_merge_base[3]), 100 # rescore needed


class _ExtendWorkerProcess(multiprocessing.Process):
    def __init__(self, anchors_to_extend: Union[str, blast_anchors_type],
                 needle_fasta_filename: str,
                 haystack_fasta_filename: str,
                 word_len: int,
                 output_queue: multiprocessing.Queue,
                 score_cutoff_slope:float,
                 score_cutoff_intersect:float
                 ):
        super().__init__()
        self.haystack_fasta = None
        self.needle_fasta = None
        self.anchors_to_extend = anchors_to_extend
        self.needle_fasta_filename = needle_fasta_filename
        self.haystack_fasta_filename = haystack_fasta_filename
        self.output_queue = output_queue
        self.logger_handler = logging.getLogger("Worker")
        self.logger_handler.debug(self.__repr__())
        self.word_len = word_len
        self.score_cutoff_slope = score_cutoff_slope
        self.score_cutoff_intersect = score_cutoff_intersect

    def extend(self, raw_anchor: blast_anchor_type) -> Optional[blast_extended_anchor_type]:
        # print(raw_anchor)
        needle_len = self.word_len
        needle_chromosome, needle_strand, needle_start = raw_anchor[0]
        needle_end = needle_start + needle_len
        haystack_chromosome, haystack_strand, haystack_start = raw_anchor[1]
        haystack_end = haystack_start + needle_len
        max_score = needle_len * 2
        max_coordinate = (needle_start, needle_end, haystack_start, haystack_end)
        needle_full_length = self.needle_fasta.get_chromosome_length(needle_chromosome)
        haystack_full_length = self.haystack_fasta.get_chromosome_length(haystack_chromosome)
        cutoff_score = self.score_cutoff_slope * needle_full_length + self.score_cutoff_intersect

        # Extend bi-direction
        left_extendable_len = min(needle_start, haystack_start)
        right_extendable_len = min(needle_full_length - needle_end, haystack_full_length - haystack_end)
        while left_extendable_len > 0 and \
                right_extendable_len > 0:
            needle_start -= 1
            haystack_start -= 1
            left_extendable_len -= 1
            needle_end += 1
            haystack_end += 1
            right_extendable_len -= 1
            needle_bytes = self.needle_fasta.get(
                needle_chromosome,
                needle_start,
                needle_end,
                needle_strand)
            haystack_bytes = self.haystack_fasta.get(
                haystack_chromosome,
                haystack_start,
                haystack_end,
                haystack_strand)
            score = simple_blast_scores(needle_bytes, haystack_bytes)
            if score > max_score:
                max_score = score
                max_coordinate = (needle_start, needle_end, haystack_start, haystack_end)
            elif score < 0:
                break
        if max_score > cutoff_score:
            return (
                (needle_chromosome, needle_strand, max_coordinate[0], max_coordinate[1] - max_coordinate[0]),
                (haystack_chromosome, haystack_strand, max_coordinate[2], max_coordinate[3] - max_coordinate[2]),
                max_score
            )
        else:
            return None

    def run(self):
        self.needle_fasta = in_memory_fasta.Fasta(self.needle_fasta_filename)
        self.haystack_fasta = in_memory_fasta.Fasta(self.haystack_fasta_filename)
        if isinstance(self.anchors_to_extend, str):
            pkl_filename = self.anchors_to_extend
            self.anchors_to_extend = pickle.load(open(pkl_filename, "rb"))
            os.remove(pkl_filename)
        for alignment in self.anchors_to_extend:
            extended_anchor = self.extend(alignment)
            if extended_anchor is not None:
                self.output_queue.put(extended_anchor)
        self.logger_handler.debug("FIN")

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
                 extend_batch_size: Optional[int] = EXTENDER_BATCH_SIZE
                 ):
        self._pool_len = pool_len
        """Internal thread number"""
        self.haystack_index = haystack_index
        self.needle_index = needle_index
        self.extend_batch_size = extend_batch_size

    def generate_raw_anchors(self) -> blast_anchors_type:
        """Generate a series of anchors"""
        haystack_available_prefixes = self.haystack_index.available_prefixes
        needle_words = self.needle_index.loaded_words
        for needle_prefix in tqdm.tqdm(desc="Anchoring...", iterable=self.needle_index.loaded_prefixes):
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
            self.haystack_index.drop_index([needle_prefix])

    def extend(self,
               raw_anchors: blast_anchors_type,
               score_cutoff_slope:float = 1.5,
               score_cutoff_intersect:float = 0.0,
               ) -> blast_extended_anchors_type:
        """"""
        this_manager = multiprocessing.Manager()
        output_queue = this_manager.Queue()
        this_batch = []
        process_pool = parallel_helper.ParallelJobQueue(
            pool_name="Extending",
            pool_size=self._pool_len,
            with_tqdm=True
        )
        worker_args = {
            "needle_fasta_filename": self.needle_index.fasta_filename,
            "haystack_fasta_filename": self.haystack_index.fasta_filename,
            "word_len": self.needle_index.word_len,
            "output_queue": output_queue,
            "score_cutoff_slope" : score_cutoff_slope,
            "score_cutoff_intersect" : score_cutoff_intersect
        }
        tmp_dir = tempfile.mkdtemp()
        n_anchor = 0
        for raw_anchor in raw_anchors:
            if len(this_batch) < self.extend_batch_size:
                this_batch.append(raw_anchor)
            else:
                pickle.dump(this_batch, open(os.path.join(tmp_dir, f"{n_anchor}.pkl"), "wb"))
                del this_batch
                gc.collect()
                this_batch = []
                n_anchor += 1
        if len(this_batch) > 0:
            pickle.dump(this_batch, open(os.path.join(tmp_dir, f"{n_anchor}.pkl"), "wb"))
            del this_batch
            gc.collect()
            n_anchor += 1
        for i in range(n_anchor):
            ewp = _ExtendWorkerProcess(anchors_to_extend=os.path.join(tmp_dir, f"{i}.pkl"), **worker_args)
            process_pool.append(ewp)

        process_pool.start()
        while not process_pool.all_finished:
            try:
                # print(ewp.exitcode, ewp.is_alive())
                yield output_queue.get(timeout=0.1)
                # print(get_item)
            except (TimeoutError, queue.Empty):
                pass
        process_pool.join()
        this_manager.shutdown()
        os.rmdir(tmp_dir)

    def extend_single_thread(self, raw_anchors: blast_anchors_type) -> blast_extended_anchors_type:
        raw_anchors = list(raw_anchors)
        this_manager = multiprocessing.Manager()
        output_queue = this_manager.Queue()
        ewp = _ExtendWorkerProcess(anchors_to_extend=raw_anchors,
                                   needle_fasta_filename=self.needle_index.fasta_filename,
                                   haystack_fasta_filename=self.haystack_index.fasta_filename,
                                   word_len=self.needle_index.word_len, output_queue=output_queue)
        ewp.start()
        with tqdm.tqdm(total=len(raw_anchors)) as pbar:
            while ewp.is_alive():
                try:
                    # print(ewp.exitcode, ewp.is_alive())
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


def _test_on_e_coli():
    needle_index = indexer.InMemorySimpleBlastIndex(word_len=5, prefix_len=1)
    needle_index.attach_fasta("test/e_coli_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/e_coli", word_len=5, prefix_len=1)
    haystack_index.attach_fasta("test/e_coli.fasta")
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index, extend_batch_size=1000)
    raw_anchors = searcher.generate_raw_anchors()
    with open("1.log", "w") as writer:
        for item in searcher.extend(raw_anchors):
            writer.write(str(item[2]) + '\n')


def _test_on_test():
    index_filename = "test/herv.pkl.xz"
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/herv.fasta")
    if os.path.exists(index_filename):
        needle_index.load(index_filename)
    else:
        needle_index.create_index()
        needle_index.save(index_filename)
    haystack_index = indexer.BlastIndex("test/test")
    haystack_index.attach_fasta("test/test.fasta")
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index)
    raw_anchors = searcher.generate_raw_anchors()
    with open("1.gtf", "w") as writer:
        for item in merge_adjacent_anchors(searcher.extend(raw_anchors)):
            writer.write(extended_anchor_to_gtf(item) + '\n')


if __name__ == "__main__":
    _test_on_test()
