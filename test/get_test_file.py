import gzip
import logging
import os.path

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
    logger_handler.info(f"Retriving {url} -> {dest_filename}")
    request_handler = requests.get(url, stream=True)
    file_size = int(request_handler.headers['content-length'])
    chunk_size = 1024

    with open(dest_filename, 'wb') as output_file:
        for chunk in tqdm(iterable=request_handler.iter_content(chunk_size=chunk_size), total=file_size // chunk_size,
                          unit='KiB'):
            if chunk:
                output_file.write(chunk)
    logger_handler.info(f"Retriving {url} -> {dest_filename} FIN")


def get_chr1_chr2_fasta():
    final_test_fa_filename = "test.fa"
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


def main():
    get_chr1_chr2_fasta()


if __name__ == '__main__':
    main()
