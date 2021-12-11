import glob
import pathlib

import setuptools
from Cython import Build
from setuptools import setup

cythonize_files = glob.glob('src/**/*.pyx', recursive=True)

PKG_NAME = "HERVfinder"

PWD = pathlib.Path(__file__).parent.resolve()

install_requires = []

with  open('requirements.txt', 'rt', encoding='utf-8') as reader:
    for line in reader:
        if not line.startswith('#'):
            install_requires.append(line.strip())
version = "1.0.0"

with  open('Readme.md', 'rt', encoding='utf-8') as reader:
    long_description = reader.read()

setup(
    name=PKG_NAME,
    version=version,
    author="YU Zhejian, XIAO Xia",
    author_email="Zhejian.19@intl.zju.edu.cn",
    description="HERVfinder -- Ultrafast BLAST-based Searcher for HERV Sequence",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="http://github.com/YU-Zhejian/HERVfinder",  # TODO
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Programming Language :: C",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: R",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.6",
    packages=setuptools.find_packages(
        where='src',
        include=['*'],
    ),
    package_dir={"": 'src'},
    install_requires=install_requires,
    ext_modules=Build.cythonize(cythonize_files),
)