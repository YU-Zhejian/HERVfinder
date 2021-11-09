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
import lzma
import pickle
from typing import Any


def dump(obj:Any, output_path:str) -> None:
    """
    Pickle the object.

    :param obj: Object that will be pickled.
    :param output_path: Output path of pickled object.
    :return: None
    """
    with lzma.open(output_path, 'wb') as writer:
        pickle.dump(obj,writer)

def load(input_path:str) -> Any:
    """
    Unpickle the object.

    :param input_path: Input path of pickled object.
    :return: The object.
    """
    with lzma.open(input_path,'rb') as reader:
        reto = pickle.load(reader)
    return reto

if __name__ == '__main__':
    import doctest
    doctest.testmod()
