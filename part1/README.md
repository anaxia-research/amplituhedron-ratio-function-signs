# Signs of Amplitudes in the Amplituhedron I: code and logs

This folder holds the computations for the paper *Signs of Amplitudes in the Amplituhedron I: The One-Loop Six-Point Ratio Function* by Deep Bhattacharjee. The paper proves that the bosonized one-loop six-point NMHV ratio function of planar N=4 super-Yang–Mills theory is strictly negative on the whole positive region. It also proves the estimate 2P/Ω ≤ −min{(u+v+w−1)/10, 10⁻³}.

## Contents

| Folder / file | Purpose |
|---|---|
| `src/sym6.py`, `src/check6.py` | Exact symbolic computation of the tree amplitude Ω, the three coefficients κ, their common denominator N, the positivity of N − n (so κ < 1), and the linear relation |
| `src/idhp.py`, `src/idcheck.py` | Check at 60 digits that the κ form equals the ratio function built from R-invariants and dilogarithms |
| `src/explore.py`, `src/lp2.py`, `src/bclosed.py`, `src/b2.py` | Floating-point versions of P, of the dual bound B and of H (used only by the independent checks) |
| `src/figdata.py` | Data for the scatter figure |
| `cert/common.py` | Interval enclosures (python-flint `arb`, 128 bits) of Li₂(1−x), μ, λ and their derivatives |
| `cert/edge.py` | Edge region e ≤ 1/100: certified slope bound |
| `cert/tails.py` | Explicit constants for the regions v ≥ e⁵ and w ≥ e²⁴ |
| `cert/small_u.py` | Two-dimensional certificate for the strip u ≤ e⁻⁷ (1,018 boxes) |
| `cert/core.py` | Three-dimensional certificate for the core region (70,171 boxes, B ≤ −0.001) |
| `cert/validate_core.py` | Independent check: 140,002 random points pushed down the same tree, no violations |
| `cert/negative_control.py` | Turns off the Δ ≥ 0 rule and shows the bisection then cannot finish, as it should |
| `cert/sanity_regions.py` | Floating-point check of every region bound on 8.3 million points |
| `cert/logs/` | Output of every run |

## Running

```
pip install python-flint sympy mpmath numpy scipy
bash run_all.sh
```

A full run takes about three minutes on one core. The certified steps take about one second each.

## Licence

MIT.
