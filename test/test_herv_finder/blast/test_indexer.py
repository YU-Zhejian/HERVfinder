import glob
import os
import pickle
import random
import tempfile
from collections import defaultdict

from herv_finder.blast.indexer import _IndexWorkerProcess, blast_index_type, _is_low_complexity, \
    _blast_index_merge

test_std_dict1 = {
    b"AAA": {("1", True, 3), ("1", True, 4)},
    b"BBB": {("1", True, 3), ("1", False, 4)}
}
test_std_dict2 = {
    b"AAA": {("1", True, 3), ("1", True, 4)},
    b"BBB": {("1", True, 3), ("1", False, 4)}
}
test_std_dict3 = {
    b"AAA": {("1", True, 3), ("1", True, 4)},
    b"CCC": {("1", True, 3), ("1", False, 4)}
}
test_std_dict4 = {
    b"AAA": {("1", True, 3), ("1", True, 4)},
    b"CCC": {("1", True, 5), ("1", False, 4)}
}

full_merged_answer = {
    b"AAA": {("1", True, 3), ("1", True, 4)},
    b"BBB": {("1", True, 3), ("1", False, 4)},
    b"CCC": {("1", True, 3), ("1", False, 4), ("1", True, 5)}
}


def test_std_merger():
    test_std_merged_dict = _blast_index_merge(test_std_dict1, test_std_dict2)
    assert test_std_dict1 == test_std_dict2 == test_std_merged_dict
    test_std_merged_dict = _blast_index_merge(test_std_merged_dict, test_std_dict3)
    test_std_merged_dict = _blast_index_merge(test_std_merged_dict, test_std_dict4)
    assert test_std_merged_dict == full_merged_answer


#
# def test_pickle_merge_with_std_dict():
#     tmp_dir = tempfile.mkdtemp()
#     pickle.dump(test_std_dict1, open(os.path.join(tmp_dir, "1.dict.pkl"), 'wb'))
#     pickle.dump(test_std_dict2, open(os.path.join(tmp_dir, "2.dict.pkl"), 'wb'))
#     pickle.dump(test_std_dict3, open(os.path.join(tmp_dir, "3.dict.pkl"), 'wb'))
#     pickle.dump(test_std_dict4, open(os.path.join(tmp_dir, "4.dict.pkl"), 'wb'))
#     new_merged_dict = {}
#     merge_pool = parallel_helper.MultiThreadingJobQueue(pool_name="Merging", with_tqdm=True)
#     merge_mutex = threading.Lock()
#
#     for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
#         merge_pool.append(_MergeWorkerThread(filename=filename, out_dict=new_merged_dict, merge_mutex=merge_mutex))
#     merge_pool.start()
#     merge_pool.join()
#     os.rmdir(tmp_dir)
#     assert new_merged_dict == full_merged_answer


def brute_force_indexer_single_strand(fasta_bytes: bytes, index_len: int) -> blast_index_type:
    tmp_dict = defaultdict(lambda: [])
    for bytes_window_start in range(len(fasta_bytes) - index_len + 1):
        # print(self.start_pos+bytes_window_start)
        bytes_window = fasta_bytes[bytes_window_start:bytes_window_start + index_len]
        if _is_low_complexity(bytes_window):
            continue

        tmp_dict[bytes_window].append(
            ('tiny', True, bytes_window_start)
        )
    return dict(tmp_dict)


def fasta_gen() -> bytes:
    chr_len = 50
    in_seq = bytearray()
    for _ in range(chr_len):
        a = random.randint(0, 3)
        in_seq.append(b'AGCT'[a])
    return bytes(in_seq)


full_length_fasta_bytes = fasta_gen()
total_len = len(full_length_fasta_bytes)
chunk_len = 10
index_len = 4


def test_same_worker_produces_same_outcome():
    tmp_dir = tempfile.mkdtemp()
    p1 = _IndexWorkerProcess(chromosome_name='tiny',
                             strand=True,
                             fasta_bytes=full_length_fasta_bytes,
                             index_len=index_len,
                             start_pos=0,
                             tmp_dir=tmp_dir)
    p2 = _IndexWorkerProcess(chromosome_name='tiny',
                             strand=True,
                             fasta_bytes=full_length_fasta_bytes,
                             index_len=index_len,
                             start_pos=0,
                             tmp_dir=tmp_dir)
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    p_outcome = []
    for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
        p_outcome.append(pickle.load(open(filename, 'rb')))
        os.remove(filename)
    merged_dict = {}
    merged_dict.update(p_outcome[0])
    merged_dict.update(p_outcome[1])
    assert p_outcome[0] == p_outcome[1]
    assert p_outcome[1] == merged_dict
    assert merged_dict == brute_force_indexer_single_strand(fasta_bytes=full_length_fasta_bytes, index_len=index_len)
    os.rmdir(tmp_dir)

# def test_merged_worker_produces_same_outcome():
#     tmp_dir = tempfile.mkdtemp()
#     pool = parallel_helper.MultiProcessingJobQueue()
#     std_new_merged_dict = {}
#     for i in range(total_len // chunk_len + 1):
#         new_worker = _IndexWorkerProcess(chromosome_name='tiny',
#                                          strand=True,
#                                          fasta_bytes=full_length_fasta_bytes[
#                                                      chunk_len * i:chunk_len * (i + 1) + index_len - 1],
#                                          index_len=index_len,
#                                          start_pos=chunk_len * i,
#                                          tmp_dir=tmp_dir)
#         pool.append(new_worker)
#     pool.start()
#     pool.join()
#
#     for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
#         _blast_index_merge(std_new_merged_dict,pickle.load(open(filename, 'rb')))
#
#     pickle_new_merged_dict = {}
#     merge_pool = parallel_helper.MultiThreadingJobQueue(pool_name="Merging", with_tqdm=True)
#     merge_mutex = threading.Lock()
#
#     for filename in glob.glob(os.path.join(tmp_dir, "*.pkl")):
#         merge_pool.append(
#             _MergeWorkerThread(filename=filename, out_dict=pickle_new_merged_dict, merge_mutex=merge_mutex))
#     merge_pool.start()
#     merge_pool.join()
#
#     os.rmdir(tmp_dir)
#     assert pickle_new_merged_dict == std_new_merged_dict == brute_force_indexer_single_strand(
#         fasta_bytes=full_length_fasta_bytes, index_len=index_len)
