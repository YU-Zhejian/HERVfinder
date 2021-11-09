"""
A general-purposed bytearray-based memory-access read-only FASTA class with GZip support.
"""
import gzip
import os.path
from typing import Optional


class Fasta:
    """
    A general-purposed bytearray-based memory-access read-only FASTA class with GZip support.

    Example:

    >>> from herv_finder.utils.in_memory_fasta import Fasta
    >>> import gzip
    >>> import os

    >>> fasta_path = "test.fa.gz"
    >>> fasta_seq = b'''>chr1 some att\\nNNNNNNNNNNNNNNNATCGTTACGTAC\\nCATATACTATATCTTAGTCTAGTCTAA\\nCGTCTTTTTCTNNNNNNNNNNNNNNNA\\nNNNNNNNNATCGTTACGTACTTCTNNN\\nCATATACTATATCTTAGTCTAGTCTAA\\nCGTCTTTTTCTNNNNNNNN\\n>chr2\\nNNNNNNNNNNNNNNNATCGTTACGTAC\\nCATATACTATATCTTAGTCTAGTCTAA\\nCGTCTTTTTCTNNNNNNNNN\\n'''
    >>> writer= gzip.open(fasta_path,'w')
    >>> _ = writer.write(fasta_seq)
    >>> writer.close()

    >>> fa = Fasta(fasta_path)
    >>> fa.get('chr1 some att', 26, 29)
    bytearray(b'CCA')
    >>> fa.get('chr2')
    bytearray(b'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN')

    >>> os.remove(fasta_path)
    """
    def __init__(self, filename:str):
        self._filename = os.path.abspath(os.path.expanduser(filename))
        """The absolute path of Fasta."""

        if not os.path.isfile(self._filename):
            raise FileNotFoundError(f"${self._filename} do not exist!")

        self._fasta_content = {}
        """
        The content of fasta. Format:  Dict[str, Array[byte]]
        """

        self._load()

    def _load(self):
        """Read the fasta file and load it into memory."""
        if self._filename.endswith(".gz") or self._filename.endswith(".GZ"):
            reader_func = gzip.open
        else:
            reader_func=open
        with reader_func(self._filename, "rb") as reader:
            chromosome_name = ""
            seq = bytearray()

            while True:
                line = reader.readline()
                if not line:
                    break
                line=line.rstrip()
                if line.startswith(b'>'):
                    if chromosome_name != "":
                        self._fasta_content[chromosome_name] = seq
                        seq = bytearray()
                    chromosome_name = str(line[1:].strip(), encoding ='UTF-8')
                else:
                    seq.extend(line)
                if chromosome_name != '':
                    self._fasta_content[chromosome_name] = seq

    def get(self, chromosome_name: str, start: Optional[int] = 0, end: Optional[int] = -1) -> bytearray:
        """Get a [) sequence for random access. If end = -1, will put all sequences till end."""
        if chromosome_name not in self._fasta_content:
            raise KeyError(f"Illegal chromosome name {chromosome_name}")
        if end == -1:
            end = len(self._fasta_content[chromosome_name])
        return self._fasta_content[chromosome_name][start:end]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
