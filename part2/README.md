# Signs of Amplitudes in the Amplituhedron II: checks and data

This folder holds the computations for the paper *Signs of Amplitudes in the Amplituhedron II: The Two-Loop Six-Point MHV Amplitude* by Deep Bhattacharjee. The paper proves that the two-loop six-point MHV amplitude of planar N=4 super-Yang–Mills theory, normalized by the BDS-like ansatz, is positive on the whole MHV positive region. The proof is analytic; the scripts here check the identities it uses, reproduce the published line expansions, and give the numerical evidence for the remainder-function conjecture.

| File | Purpose |
|---|---|
| `e2lib.py` | The two-loop amplitude in the cube coordinates (a, b, c): the five-term decomposition of E^(2) and R6^(2) |
| `r62.py` | The Goncharov–Spradlin–Vergu–Volovich formula in the cross ratios u, v, w |
| `bij.py` | Symbolic proof of the rational parametrization of the positive region by the unit cube, and of the product identity behind J = −δ/2 |
| `identities.py` | 40-digit checks: decomposition against the GSVV formula, J = −δ/2, the derivative formula for L4, the reflection identity |
| `linecheck.py` | Reproduces the two-loop line coefficients of Dixon–von Hippel–McLeod–Trnka at u = 0.3 |
| `evidence.py` | 40,000-point sample of the cube: R6^(2) > 0 at every point (conjecture), smallest values of E^(2) and Y |
| `figdata.py` | Data for the figures (`figdata/`) |
| `logs/` | Output of every run |

## Running

```
pip install mpmath sympy
bash run_all.sh
```

A full run takes about three minutes on four cores.
