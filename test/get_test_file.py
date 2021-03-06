import gzip
import logging
import os.path
from typing import BinaryIO

import requests
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger_handler = logging.getLogger()


def download_file(dest_filename: str, url: str):
    """
    Download some file from url to dest_filename

    This code uses codes in https://www.simplifiedpython.net/python-download-file/
    """
    global logger_handler
    logger_handler.info(f"Retrieving {url} -> {dest_filename}")
    request_handler = requests.get(url, stream=True)
    file_size = int(request_handler.headers['content-length'])
    chunk_size = 1024

    with open(dest_filename, 'wb') as output_file:
        for chunk in tqdm(iterable=request_handler.iter_content(chunk_size=chunk_size), total=file_size // chunk_size,
                          unit='KiB'):
            if chunk:
                output_file.write(chunk)
    logger_handler.info(f"Retrieving {url} -> {dest_filename} FIN")


def download_stream(writer: BinaryIO, url: str):
    """
    Download some file from url to dest_filename

    This code uses codes in https://www.simplifiedpython.net/python-download-file/
    """
    global logger_handler
    logger_handler.info(f"Retrieving {url}")
    request_handler = requests.get(url, stream=True)
    file_size = int(request_handler.headers['content-length'])
    chunk_size = 1024

    for chunk in tqdm(iterable=request_handler.iter_content(chunk_size=chunk_size), total=file_size // chunk_size,
                      unit='KiB'):
        if chunk:
            writer.write(chunk)
    logger_handler.info(f"Retrieving {url} FIN")


def get_chr21_chr22_fasta():
    final_test_fa_filename = "test.fasta"
    if os.path.exists(final_test_fa_filename):
        return
    with open(final_test_fa_filename, 'wb') as writer:
        for chromosome_name in ("chr21", "chr22"):
            gzipped_filename = f"{chromosome_name}.fa.gz"
            if not os.path.exists(gzipped_filename):
                download_file(gzipped_filename,
                              f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{gzipped_filename}")
            with gzip.open(gzipped_filename, 'rb') as reader:
                while True:
                    contents_to_write = reader.read(1024)
                    if not contents_to_write:
                        break
                    writer.write(contents_to_write)


def get_e_coli_fasta():
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz'
    final_test_fa_filename = "e_coli.fasta"
    if os.path.exists(final_test_fa_filename):
        return
    with open(final_test_fa_filename, 'wb') as writer:
        gzipped_filename = "e_coli.fasta.gz"
        if not os.path.exists(gzipped_filename):
            download_file(gzipped_filename, url)
            with gzip.open(gzipped_filename, 'rb') as reader:
                while True:
                    contents_to_write = reader.read(1024)
                    if not contents_to_write:
                        break
                    writer.write(contents_to_write)


def get_herv_seq():
    url = 'https://www.dfam.org/api/families?format=fasta&name_accession=LTR&clade=9606&clade_relatives=descendants&download=true'
    final_test_fa_filename = "herv.fasta"
    if not os.path.exists(final_test_fa_filename):
        download_file(final_test_fa_filename, url)


def get_herv_tsv():
    urls = [
        'https://www.dfam.org/api/families/DF0000472/assemblies/hg38/annotations?nrph=false&download=true',
        'https://www.dfam.org/api/families/DF0000558/assemblies/hg38/annotations?nrph=false&download=true',
        'https://www.dfam.org/api/families/DF0000471/assemblies/hg38/annotations?nrph=false&download=true'
    ]
    final_test_fa_filename = "herv.tsv"
    if not os.path.exists(final_test_fa_filename):
        with open(final_test_fa_filename, 'wb') as writer:
            for url in urls:
                download_stream(writer, url)


def main():
    os.chdir(os.path.dirname(__file__))
    get_chr21_chr22_fasta()
    get_e_coli_fasta()
    get_herv_seq()
    get_herv_tsv()


if __name__ == '__main__':
    main()
