"""
A general-purposed bytes-based memory-access read-only FASTA class with GZip support.
"""
import gzip
import logging
import os.path
from typing import Optional, List

_FASTA_COMP_TRANS = bytes.maketrans(b'ATCGatcgNnXx', b'TAGCtagcNnXx')
"""The fasta complementary translator dictionary"""

_FASTA_REMOVE_LOWCASE_TRANS = bytes.maketrans(b'atcgnx', b'ATCGNX')
"""The fasta low-to-high translator dictionary"""

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()


class Fasta:
    """
    A general-purposed bytes-based memory-access read-only FASTA class with GZip support.

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
    b'CCA'
    >>> fa.get('chr2')
    b'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN'

    >>> os.remove(fasta_path)
    """

    def __init__(self, filename: str):
        logger_handler.info(f"Creating FASTA from {filename}...")
        self._filename = os.path.abspath(os.path.expanduser(filename))
        """The absolute path of Fasta."""

        if not os.path.isfile(self._filename):
            raise FileNotFoundError(f"${self._filename} do not exist!")

        self._fasta_content = {}
        """
        The content of fasta. Format:  Dict[str, bytes]
        """

        self._load()
        logger_handler.info(f"FASTA load complete")

    def _load(self):
        """Read the fasta file and load it into memory."""
        if self._filename.endswith(".gz") or self._filename.endswith(".GZ"):
            reader_func = gzip.open
        else:
            reader_func = open
        with reader_func(self._filename, "rb") as reader:
            chromosome_name = ""
            seq = bytearray()
            while True:
                line = reader.readline()
                if not line:
                    break
                line = line.rstrip()
                if line.startswith(b'>'):
                    if chromosome_name != "":
                        logger_handler.debug(f"Chromosome {chromosome_name} FIN")
                        self._fasta_content[chromosome_name] = bytes(seq).translate(_FASTA_REMOVE_LOWCASE_TRANS)

                        seq = bytearray()
                    chromosome_name = str(line[1:].strip(), encoding='UTF-8')
                    logger_handler.debug(f"New chromosome {chromosome_name}")
                else:
                    seq.extend(line)
            if chromosome_name != '':
                logger_handler.debug(f"Chromosome {chromosome_name} FIN")
                self._fasta_content[chromosome_name] = bytes(seq).translate(_FASTA_REMOVE_LOWCASE_TRANS)

    def get(self, chromosome_name: str, start: Optional[int] = 0, end: Optional[int] = -1,
            strand: bool = True) -> bytes:
        """Get a 0-based [) sequence for random access. If end = -1, will put all sequences till end."""

        if chromosome_name not in self._fasta_content:
            raise KeyError(f"Illegal chromosome name {chromosome_name}")
        if end == -1:
            end = len(self._fasta_content[chromosome_name])
        if strand:
            return self._fasta_content[chromosome_name][start:end]
        else:
            return get_reversed_complementary(self._fasta_content[chromosome_name][start:end])

    @property
    def chromosomes(self) -> List[str]:
        return list(self._fasta_content.keys())


def get_reversed_complementary(fasta_bytes: bytes) -> bytes:
    return fasta_bytes.translate(_FASTA_COMP_TRANS)[::-1]


if __name__ == "__main__":
    import doctest

    doctest.testmod()
