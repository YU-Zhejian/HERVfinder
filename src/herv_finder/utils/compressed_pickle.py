"""
Pickle with LZMA compression.

Example:

>>> import os
>>> from herv_finder.utils import compressed_pickle
>>> some_object = "Some text"
>>> pickle_path = "pickled.pkl"
>>> compressed_pickle.dump(some_object, pickle_path)
>>> unpickled_some_object = compressed_pickle.load(pickle_path)
>>> assert some_object == unpickled_some_object
>>> os.remove(pickle_path)
"""
import io
import lzma
import os
import pickle
from typing import Any

import tqdm


def dump(obj: Any, output_path: str) -> None:
    """
    Pickle the object.

    :param obj: Object that will be pickled.
    :param output_path: Output path of pickled object.
    :return: None
    """
    with lzma.open(output_path, 'wb', preset=lzma.PRESET_EXTREME) as writer:
        pickle.dump(obj, writer)


def load(input_path: str) -> Any:
    """
    Unpickle the object.

    :param input_path: Input path of pickled object.
    :return: The object.
    """
    with lzma.open(input_path, 'rb') as reader:
        reto = pickle.load(reader)
    return reto


class _UnpickleLZMAWithTQDM:
    """
    The underlying class with standard functions defined.
    """

    def __init__(self, path, **kwargs):
        super().__init__()
        self.fd = lzma.open(path, 'rb')
        curr_pos = self.fd.tell()
        self.fd.seek(offset=0, whence=io.SEEK_END)
        total_size = self.fd.tell()
        self.fd.seek(curr_pos)
        self.tqdm = tqdm.tqdm(desc=f"Unpickling {path}", total=total_size, **kwargs)

    def readline(self):
        update_bytes = self.fd.readline()
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def read(self, size: int = -1):
        update_bytes = self.fd.read(size)
        self.tqdm.update(len(update_bytes))
        return update_bytes

    def __enter__(self):
        self.tqdm.__enter__()
        return self

    def __exit__(self, *args, **kwargs):
        return self.tqdm.__exit__(*args, **kwargs)


# noinspection all
def load_with_tqdm(filename: str):
    """
    Unpickle a file with tqdm.

    :return: Picked object.
    """
    with _UnpickleLZMAWithTQDM(filename) as pbfd:
        up = pickle.Unpickler(pbfd)
        obj = up.load()
    return obj


if __name__ == '__main__':
    import doctest

    doctest.testmod()
