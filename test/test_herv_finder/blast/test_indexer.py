from herv_finder.blast.indexer import _blast_index_merge

test_std_dict1 = {
    b"BBB": [("1", True, 3), ("1", False, 4)],
    b"CCC": [("1", True, 6), ("1", False, 5)]
}
test_std_dict2 = {
    b"AAA": [("1", True, 3), ("1", True, 4)],
    b"CCC": [("1", True, 3), ("1", False, 4)]
}


def test_merge():
    _blast_index_merge(test_std_dict1, test_std_dict2)
    assert test_std_dict1 == {b'BBB': [('1', True, 3), ('1', False, 4)],
                              b'CCC': [('1', True, 6), ('1', False, 5), ("1", True, 3), ("1", False, 4)],
                              b'AAA': [('1', True, 3), ('1', True, 4)]}
