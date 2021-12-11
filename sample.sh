#!/usr/bin/env bash
# A sample Shell script
set -ue
cd "$(dirname "$(readlink -f "${0}")")" || exit

PYTHONPATH="${PYTHONPATH:-}:src" python -m herv_finder index -R test/test.fasta -I test/test
PYTHONPATH="${PYTHONPATH:-}:src" python -m herv_finder search -R test/test.fasta -I test/test -H test/herv.fasta -O ./
