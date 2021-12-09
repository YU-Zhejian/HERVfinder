import gc
import itertools
import logging
import multiprocessing
import os
import pickle
import queue
import tempfile
from typing import Optional, Iterable

import numpy as np
import tqdm

from herv_finder.blast import indexer, DEFAULT_SCORE_CUTOFF_SLOPE, DEFAULT_SCORE_CUTOFF_INTERSECT, \
    DEFAULT_EXTENDER_BATCH_SIZE, BlastAnchor
from herv_finder.utils import in_memory_fasta, parallel_helper


def _hamming_distance_blast_score(b1: bytes, b2: bytes) -> int:
    """A very simple blast score by hamming distance"""
    assert len(b1) == len(b2)
    score = 2 * len(b1)
    for i_b, j_b in zip(b1, b2):
        if i_b != j_b:
            score -= 1
    return score


def _editing_distance_score(b1: bytes, b2: bytes) -> int:
    """A complicated blast score by editing distance"""
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


SCORING_FUNCTIONS = {
    "Hamming": _hamming_distance_blast_score,
    "Edit Distance": _editing_distance_score
}


def rescore(needle_fasta: in_memory_fasta.Fasta,
            haystack_fasta: in_memory_fasta.Fasta,
            anchor: BlastAnchor,
            scoring_function: str = "Edit Distance"
            ) -> int:
    needle_bytes = needle_fasta.get(
        anchor.needle_chromosome,
        anchor.needle_start,
        anchor.needle_end,
        anchor.strand)
    haystack_bytes = haystack_fasta.get(
        anchor.haystack_chromosome,
        anchor.haystack_start,
        anchor.haystack_end,
        anchor.strand)
    return SCORING_FUNCTIONS[scoring_function](needle_bytes, haystack_bytes)


class _ExtendWorkerProcess(multiprocessing.Process):
    anchors_to_extend: Iterable[BlastAnchor]
    haystack_fasta: in_memory_fasta.Fasta
    needle_fasta: in_memory_fasta.Fasta

    def __init__(self, anchors_to_extend_filename: str,
                 word_len: int,
                 needle_fasta_filename: str,
                 haystack_fasta_filename: str,
                 output_queue: multiprocessing.Queue,
                 score_cutoff_slope: float,
                 score_cutoff_intersect: float
                 ):
        super().__init__()
        self.anchors_to_extend_filename = anchors_to_extend_filename
        self.needle_fasta_filename = needle_fasta_filename
        self.haystack_fasta_filename = haystack_fasta_filename
        self.output_queue = output_queue
        self.logger_handler = logging.getLogger("Worker")
        self.logger_handler.debug(self.__repr__())
        self.word_len = word_len
        self.score_cutoff_slope = score_cutoff_slope
        self.score_cutoff_intersect = score_cutoff_intersect

    def extend(self, anchor: BlastAnchor) -> Optional[BlastAnchor]:
        max_score = anchor.needle_len + anchor.haystack_len
        needle_full_length = self.needle_fasta.get_chromosome_length(anchor.needle_chromosome)
        haystack_full_length = self.haystack_fasta.get_chromosome_length(anchor.haystack_chromosome)
        cutoff_score = self.score_cutoff_slope * needle_full_length + self.score_cutoff_intersect

        # Extend bi-direction
        max_start = (anchor.needle_start, anchor.haystack_start)
        max_end = (anchor.needle_end, anchor.haystack_end)
        while anchor.needle_start > 0 and anchor.haystack_start > 0:
            anchor.needle_start -= 1
            anchor.haystack_start -= 1
            score = rescore(self.needle_fasta, self.haystack_fasta, anchor)
            if score > max_score:
                max_score = score
                max_start = (anchor.needle_start, anchor.haystack_start)
            elif score < 0:
                break
        anchor.needle_start, anchor.haystack_start = max_start
        while anchor.needle_end < needle_full_length and anchor.haystack_end < haystack_full_length:
            anchor.needle_end += 1
            anchor.haystack_end += 1
            score = rescore(self.needle_fasta, self.haystack_fasta, anchor)
            if score > max_score:
                max_score = score
                max_end = (anchor.needle_end, anchor.haystack_end)
            elif score < 0:
                break
        anchor.needle_end, anchor.haystack_end = max_end
        if max_score > cutoff_score:
            anchor.score = max_score
            return anchor
        else:
            return None

    def run(self):
        self.needle_fasta = in_memory_fasta.Fasta(self.needle_fasta_filename)
        self.haystack_fasta = in_memory_fasta.Fasta(self.haystack_fasta_filename)
        self.anchors_to_extend = pickle.load(open(self.anchors_to_extend_filename, "rb"))
        os.remove(self.anchors_to_extend_filename)
        for anchor in self.anchors_to_extend:
            extended_anchor = self.extend(anchor)
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
                 extend_batch_size: Optional[int] = DEFAULT_EXTENDER_BATCH_SIZE
                 ):
        self._pool_len = pool_len
        """Internal thread number"""
        self.haystack_index = haystack_index
        self.needle_index = needle_index
        self.extend_batch_size = extend_batch_size

    def _generate_raw_anchors(self) -> Iterable[BlastAnchor]:
        """Generate a series of anchors"""
        haystack_available_prefixes = self.haystack_index.available_prefixes
        needle_words = self.needle_index.loaded_words
        with open("raw_anchor.txt", "w") as writer:
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
                        raw_anchor = BlastAnchor.from_location(needle_location, haystack_location,
                                                               self.needle_index.word_len)
                        writer.write(str(raw_anchor) + "\n")
                        yield raw_anchor
            self.haystack_index.drop_index([needle_prefix])

    def extend(self,
               raw_anchors: Iterable[BlastAnchor],
               score_cutoff_slope: float = DEFAULT_SCORE_CUTOFF_SLOPE,
               score_cutoff_intersect: float = DEFAULT_SCORE_CUTOFF_INTERSECT,
               ) -> Iterable[BlastAnchor]:
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
            "score_cutoff_slope": score_cutoff_slope,
            "score_cutoff_intersect": score_cutoff_intersect
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
            ewp = _ExtendWorkerProcess(anchors_to_extend_filename=os.path.join(tmp_dir, f"{i}.pkl"), **worker_args)
            process_pool.append(ewp)
        process_pool.start()
        with open("unfiltered_extended_anchor.txt", "w") as writer:
            while not process_pool.all_finished:
                try:
                    # print(ewp.exitcode, ewp.is_alive())
                    extended_anchor = output_queue.get(timeout=0.1)
                    writer.write(str(extended_anchor) + "\n")
                    yield extended_anchor
                    # print(get_item)
                except (TimeoutError, queue.Empty):
                    pass
        process_pool.join()
        this_manager.shutdown()
        os.rmdir(tmp_dir)


def _test_tiny():
    needle_index = indexer.InMemorySimpleBlastIndex()
    needle_index.attach_fasta("test/tiny_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/tiny")
    haystack_index.attach_fasta("test/tiny.fasta")
    haystack_index.create_index()
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index)
    for anchor in searcher._generate_raw_anchors():
        print(anchor)


def _test_on_e_coli():
    needle_index = indexer.InMemorySimpleBlastIndex(word_len=5, prefix_len=1)
    needle_index.attach_fasta("test/e_coli_needle.fasta")
    needle_index.create_index()
    haystack_index = indexer.BlastIndex("test/e_coli", word_len=5, prefix_len=1)
    haystack_index.attach_fasta("test/e_coli.fasta")
    searcher = BlastIndexSearcher(needle_index=needle_index, haystack_index=haystack_index, extend_batch_size=1000)
    raw_anchors = searcher._generate_raw_anchors()
    with open("1.log", "w") as writer:
        for item in searcher.extend(raw_anchors):
            writer.write(str(item) + '\n')


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
    raw_anchors = searcher._generate_raw_anchors()
    with open("1.gtf", "w") as writer:
        for item in searcher.extend(raw_anchors):
            writer.write(item.to_gtf() + '\n')


if __name__ == "__main__":
    _test_on_test()
