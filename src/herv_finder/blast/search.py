import gc
import itertools
import logging
import multiprocessing
import os
import pickle
import queue
import tempfile
from collections import defaultdict
from typing import Optional, Iterable

import tqdm

from herv_finder.blast import indexer, DEFAULT_SCORE_CUTOFF_SLOPE, DEFAULT_SCORE_CUTOFF_INTERSECT, \
    DEFAULT_EXTENDER_BATCH_SIZE, BlastAnchor
from herv_finder.blast._cblast import _editing_distance_blast_score, _hamming_distance_blast_score
from herv_finder.utils import in_memory_fasta, parallel_helper

SCORING_FUNCTIONS = {
    "Hamming": _hamming_distance_blast_score,
    "Edit Distance": _editing_distance_blast_score
}


def rescore(needle_fasta: in_memory_fasta.Fasta,
            haystack_fasta: in_memory_fasta.Fasta,
            anchor: BlastAnchor,
            scoring_function=_editing_distance_blast_score
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
    return scoring_function(needle_bytes, haystack_bytes)


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
                 score_cutoff_intersect: float,
                 output_basename: str,
                 io_mutex: multiprocessing.Lock
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
        self.io_mutex = io_mutex
        self.output_basename = output_basename

    def pre_extend_merge(self):
        """"""
        pass

    def extend(self, anchor: BlastAnchor) -> Optional[BlastAnchor]:
        # print(anchor)
        max_score = 0
        needle_full_length = self.needle_fasta.get_chromosome_length(anchor.needle_chromosome)
        haystack_full_length = self.haystack_fasta.get_chromosome_length(anchor.haystack_chromosome)
        cutoff_score = self.score_cutoff_slope * needle_full_length + self.score_cutoff_intersect

        # Extend bi-direction
        max_start = (anchor.needle_start, anchor.haystack_start)
        max_end = (anchor.needle_end, anchor.haystack_end)
        while anchor.needle_start > 0 and anchor.haystack_start > 0:
            anchor.needle_start -= 1
            anchor.haystack_start -= 1
            score = rescore(self.needle_fasta, self.haystack_fasta, anchor, _hamming_distance_blast_score)
            if score >= max_score:
                max_score = score
                max_start = (anchor.needle_start, anchor.haystack_start)
            elif score < 0:
                break
        anchor.needle_start, anchor.haystack_start = max_start
        max_score = 0
        while anchor.needle_end < needle_full_length and anchor.haystack_end < haystack_full_length:
            anchor.needle_end += 1
            anchor.haystack_end += 1
            score = rescore(self.needle_fasta, self.haystack_fasta, anchor, _hamming_distance_blast_score)
            if score >= max_score:
                max_score = score
                max_end = (anchor.needle_end, anchor.haystack_end)
            elif score < 0:
                break
        anchor.needle_end, anchor.haystack_end = max_end
        anchor.score = max_score
        with self.io_mutex:
            with open(self.output_basename + "unfiltered_extended_anchor.txt", "a") as txt_writer:
                with open(self.output_basename + "unfiltered_extended_anchor.gtf", "a") as gtf_writer:
                    txt_writer.write(str(anchor) + "\n")
                    txt_writer.flush()
                    gtf_writer.write(anchor.to_gtf() + "\n")
                    gtf_writer.flush()
        if max_score > cutoff_score:
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
                 extend_batch_size: Optional[int] = DEFAULT_EXTENDER_BATCH_SIZE,
                 output_basename: Optional[str] = ""
                 ):
        self._pool_len = pool_len
        """Internal thread number"""
        self.haystack_index = haystack_index
        self.needle_index = needle_index
        self.extend_batch_size = extend_batch_size
        self.needle_fasta = in_memory_fasta.Fasta(self.needle_index.fasta_filename)
        self.haystack_fasta = in_memory_fasta.Fasta(self.haystack_index.fasta_filename)
        self.output_basename = output_basename

    def generate_raw_anchors(self) -> Iterable[BlastAnchor]:
        """Generate a series of anchors"""
        haystack_available_prefixes = self.haystack_index.available_prefixes
        needle_words = self.needle_index.loaded_words
        with open(self.output_basename + "raw_anchor.txt", "w") as txt_writer:
            with  open(self.output_basename + "raw_anchor.gtf", "w") as gtf_writer:
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
                            txt_writer.write(str(raw_anchor) + "\n")
                            gtf_writer.write(raw_anchor.to_gtf() + "\n")
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
        io_mutex = this_manager.Lock()
        worker_args = {
            "needle_fasta_filename": self.needle_index.fasta_filename,
            "haystack_fasta_filename": self.haystack_index.fasta_filename,
            "word_len": self.needle_index.word_len,
            "output_queue": output_queue,
            "score_cutoff_slope": score_cutoff_slope,
            "score_cutoff_intersect": score_cutoff_intersect,
            "output_basename": self.output_basename,
            "io_mutex": io_mutex
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
        try:
            os.remove(self.output_basename + "unfiltered_extended_anchor.txt")
            os.remove(self.output_basename + "unfiltered_extended_anchor.gtf")
        except FileNotFoundError:
            pass
        for i in range(n_anchor):
            ewp = _ExtendWorkerProcess(anchors_to_extend_filename=os.path.join(tmp_dir, f"{i}.pkl"), **worker_args)
            process_pool.append(ewp)
        process_pool.start()
        with open(self.output_basename + "filtered_extended_anchor.txt", "w") as txt_writer:
            with open(self.output_basename + "filtered_extended_anchor.gtf", "w") as gtf_writer:
                while not process_pool.all_finished:
                    try:
                        # print(ewp.exitcode, ewp.is_alive())
                        extended_anchor = output_queue.get(timeout=0.1)
                        txt_writer.write(str(extended_anchor) + "\n")
                        gtf_writer.write(extended_anchor.to_gtf() + "\n")
                        yield extended_anchor
                        # print(get_item)
                    except (TimeoutError, queue.Empty):
                        pass
        process_pool.join()
        this_manager.shutdown()
        os.rmdir(tmp_dir)

    def merge_overlapping_anchors(self, raw_anchors: Iterable[BlastAnchor]) -> Iterable[BlastAnchor]:
        anchor_dict = defaultdict(lambda: [])
        for anchor in raw_anchors:
            anchor_dict[(anchor.haystack_chromosome, anchor.needle_chromosome, anchor.strand)].append(anchor)
        with open(self.output_basename + "merged_anchor.txt", "w") as txt_writer:
            with open(self.output_basename + "merged_anchor.gtf", "w") as gtf_writer:
                for v in anchor_dict.values():
                    v = list(sorted(v, key=lambda x: x.needle_start))
                    v = list(sorted(v, key=lambda x: x.haystack_start))
                    current_anchor = v.pop()
                    while len(v) > 0:
                        this_anchor = v.pop()

                        if (
                                current_anchor.haystack_start <= this_anchor.haystack_start <= current_anchor.haystack_end and
                                current_anchor.needle_start <= this_anchor.needle_start <= current_anchor.needle_start
                        ) or (
                                current_anchor.haystack_start <= this_anchor.haystack_end <= current_anchor.haystack_end and
                                current_anchor.needle_start <= this_anchor.needle_end <= current_anchor.needle_start
                        ) or (
                                this_anchor.haystack_start <= current_anchor.haystack_start and
                                current_anchor.haystack_end <= this_anchor.haystack_end and
                                this_anchor.needle_start <= current_anchor.needle_start and
                                current_anchor.needle_end <= this_anchor.needle_end
                        ):
                            current_anchor.haystack_start = min(current_anchor.haystack_start,
                                                                this_anchor.haystack_start)
                            current_anchor.haystack_end = max(current_anchor.haystack_end, this_anchor.haystack_end)
                            current_anchor.needle_start = min(current_anchor.needle_start, this_anchor.needle_start)
                            current_anchor.needle_end = max(current_anchor.needle_end, this_anchor.needle_end)

                        else:
                            current_anchor.score = rescore(self.needle_fasta, self.haystack_fasta,
                                                           current_anchor,
                                                           scoring_function=_editing_distance_blast_score)
                            yield current_anchor
                            txt_writer.write(str(current_anchor) + "\n")
                            gtf_writer.write(current_anchor.to_gtf() + "\n")
                            current_anchor = this_anchor
                    current_anchor.score = rescore(self.needle_fasta, self.haystack_fasta,
                                                   current_anchor, scoring_function=_editing_distance_blast_score)
                    yield current_anchor
                    txt_writer.write(str(current_anchor) + "\n")
                    gtf_writer.write(current_anchor.to_gtf() + "\n")
