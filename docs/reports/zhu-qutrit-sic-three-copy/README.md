# Exact qutrit SIC three-copy fidelity checks

This checks the proposed family calculation for Huangjun Zhu, *Quantum Measurements in the Light of Quantum State Estimation*, PRX Quantum 3, 030306 (2022), DOI [10.1103/PRXQuantum.3.030306](https://doi.org/10.1103/PRXQuantum.3.030306), arXiv:2111.02841v2, Conjecture 1. Source quotations in the preregistration are attributed to Zhu (journal CC BY 4.0). This directory contains independently written checking code; no source PDF or source TeX is redistributed.

Requirements: Python 3.12 and SymPy 1.14.0 (`python3 -m pip install sympy==1.14.0`). Run from any directory:

```sh
python3 /path/to/check.py
```

`check.py` uses exact rational Laurent polynomials modulo `w^2+w+1`, with conjugation `t -> 1/t`, `w -> w^2`. It verifies the SIC equations, all 84 unordered distinct triples, the source operator against an independent 24-permutation contraction, characteristic polynomials for every repeated-outcome type, phase multiplicities, spectral identities, harmonic cancellation and normalization. It independently reproduces the source two-copy qutrit SIC fidelity `17/30`. The three-copy formula and endpoint comparison are exact identities. The final line is `ALL_CHECKS_OK`.

No floating-point values support the checks. `output.txt` is the complete successful output. Positivity and branch choices require the analytic argument in the preregistration. This is neither a Lean proof nor verification of the published classification of arbitrary qutrit SICs. In particular, the finite enumeration of triples in a symbolic family does not establish that all qutrit SICs belong to that family.

The script has no path-dependent inputs and was also run from a different working directory and a directory containing spaces, using a shell without login configuration. Verification platform: macOS, Python 3.12, SymPy 1.14.0. Other platforms have not been tested.
