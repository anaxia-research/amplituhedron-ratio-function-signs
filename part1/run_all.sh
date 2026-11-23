#!/usr/bin/env bash
# Reproduces every computation of Part I. Requires Python 3 with
# python-flint, sympy, mpmath, numpy and scipy.
set -euo pipefail
cd "$(dirname "$0")"

echo "== symbolic identities (R-invariants, kappa numerators, linear relation)"
( cd src && python3 sym6.py && python3 check6.py && python3 idhp.py )

echo "== certified bounds (interval arithmetic, logs in cert/logs/)"
cd cert
python3 edge.py
python3 tails.py
python3 small_u.py
python3 core.py

echo "== independent checks"
python3 validate_core.py | tee logs/validate_core.log
python3 negative_control.py
python3 sanity_regions.py | tee logs/sanity_regions.log
