


This project is designed and tested in Debian GNU/Linux 11, but it may also run in other systems.

Some core functions is accelerated using Cython. This needs you to have a working C Compiler, like GCC, Intel oneAPI C Compiler, Clang or MSVC.

## Building for Own Purpose

Make sure the virtual environments are properly configured.

Use the following command to build the Cython extension:

```shell
cythonize -i src/herv_finder/blast/_cblast.pyx
```

Then you may run this project using:

```shell
PYTHONPATH="${PYTHONPATH:-}:src" python -m herv_finder [ARGUMENTS]
```

to run this project. A sample shell script is attached in `sample.sh`.

## Building for Distribution

Make sure the virtual environments are properly configured.

If you have GNU Make installed on your machine (BSD Make will not work), type:

```shell
make dist
```

or manually invoke:

```shell
python setup.py sdist bdist_wheel
```

You may distribute the sdist (A GNU TAR archive in `dist/`, by default) to others and let them install this software using:

```shell
python -m pip install ${PATH_TO_SDIST}
```

The binary wheel may be audited by PYPA `manylinux` project to make it ABI-independent. Otherwise, it will be ABI-dependent to your CPython and GLibc.
