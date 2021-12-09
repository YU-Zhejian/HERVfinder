import argparse
import os
import sys
from typing import List

from herv_finder import blast
from herv_finder.blast import indexer, search

PROG = "HERVfinder"

BANNER = """
    =========================================================
    =                       HERVfinder                      =
    =========================================================
    """


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=PROG)
    parser.add_argument('action', choices=['index', 'search'], help='Action to do')
    parser.add_argument('-H', '--herv_fasta', type=str, help="[search] Fasta for HERV consensus sequence",
                        required=False)
    parser.add_argument('-R', '--reference_fasta', type=str,
                        help="Fasta for (a subset of) human reference genome sequence", required=True)
    parser.add_argument('-I', '--index', type=str, help="Basename of HERVfinder Blast index", required=True)
    parser.add_argument('-P', '--process', type=int, help="Number of process to use", required=False,
                        default=os.cpu_count())
    parser.add_argument('--prefix_len', type=int, help="Prefix length of splitted index", required=False,
                        default=blast.DEFAULT_PREFIX_LEN)
    parser.add_argument('--word_len', type=int, help="Word length of index", required=False,
                        default=blast.DEFAULT_WORD_LEN)
    parser.add_argument('--chunk_len', type=int, help="[index] Chunk length of each process", required=False,
                        default=blast.DEFAULT_CHUNK_LEN)
    parser.add_argument('-O', '--output', type=str, help="[search] Basename of the output", required=False)
    parser.add_argument('--score_cutoff_slope', type=float, help="[search] Scope of post-alignment score cutoff",
                        required=False, default=blast.DEFAULT_SCORE_CUTOFF_SLOPE)
    parser.add_argument('--score_cutoff_intersect', type=float,
                        help="[search] Intersect of post-alignment score cutoff",
                        required=False, default=blast.DEFAULT_SCORE_CUTOFF_INTERSECT)
    parser.add_argument('--extend_batch_size', type=int, help="[search] Number of anchors to submit to each process",
                        required=False, default=blast.DEFAULT_EXTENDER_BATCH_SIZE)
    return parser.parse_args(args)


def _search(
        herv_fasta: str,
        reference_fasta: str,
        index: str,
        output: str,
        score_cutoff_slope: float,
        score_cutoff_intersect: float,
        pool_len: int,
        extend_batch_size: int,
        prefix_len: int,
        word_len: int
):
    needle_index = indexer.InMemorySimpleBlastIndex(word_len=word_len, prefix_len=prefix_len)
    needle_index.attach_fasta(herv_fasta)
    haystack_index = indexer.BlastIndex(index, word_len=word_len, prefix_len=prefix_len)
    haystack_index.attach_fasta(reference_fasta)
    searcher = search.BlastIndexSearcher(
        needle_index=needle_index,
        haystack_index=haystack_index,
        output_basename=output,
        pool_len=pool_len,
        extend_batch_size=extend_batch_size
    )
    _ = list(searcher.merge_overlapping_anchors(searcher.extend(
        searcher.generate_raw_anchors(),
        score_cutoff_slope=score_cutoff_slope,
        score_cutoff_intersect=score_cutoff_intersect
    )))


def _index(
        reference_fasta: str,
        index: str,
        pool_len: int,
        prefix_len: int,
        word_len: int,
        chunk_len: int
):
    bi = indexer.BlastIndex(
        basename=index,
        word_len=word_len,
        chunk_len=chunk_len,
        pool_len=pool_len,
        prefix_len=prefix_len
    )
    bi.attach_fasta(reference_fasta)
    bi.create_index()
    bi.validate_index()


def main(args: List[str]) -> int:
    args = _parse_args(args)
    if args.action == 'index':
        _index(
            reference_fasta=args.reference_fasta,
            index=args.index,
            pool_len=args.pool_len,
            prefix_len=args.prefix_len,
            word_len=args.word_len,
            chunk_len=args.chunk_len
        )
    elif args.action == 'search':
        _search(
            herv_fasta=args.herv_fasta,
            reference_fasta=args.reference_fasta,
            index=args.index,
            output=args.output,
            score_cutoff_slope=args.score_cutoff_slope,
            score_cutoff_intersect=args.score_cutoff_intersect,
            pool_len=args.pool_len,
            extend_batch_size=args.extend_batch_size,
            prefix_len=args.prefix_len,
            word_len=args.word_len,
        )
    else:
        print('Undefined behaviour -- Your argparse module may not be working.')
    return 0


if __name__ == '__main__':
    print(BANNER)
    print("RECV:" + " ".join(sys.argv))
    sys.exit(main(sys.argv[1:]))
