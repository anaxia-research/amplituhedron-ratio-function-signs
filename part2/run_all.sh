#!/usr/bin/env bash
# Reproduces every check of Part II. Requires Python 3 with mpmath and sympy.
set -euo pipefail
cd "$(dirname "$0")"
python3 bij.py        | tee logs/bij.log
python3 identities.py
python3 linecheck.py  | tee logs/linecheck.log
python3 evidence.py
python3 figdata.py    > logs/figdata.log
