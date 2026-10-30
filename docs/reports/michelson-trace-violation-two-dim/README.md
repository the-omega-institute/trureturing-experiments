# Two-dimensional spectral-mean trace violation

The target is exactly item 6 of Conjecture 11.1 in Abed–Nikolaeva–Novikov,
*How characterisations of traces and central elements fail: sharp violation
constants via Michelson contrast and noncentrality*, arXiv:2610.34388v1.
The proposed answer is V₂(δ) = (1 − δ²)^(−1/2) − 1 for 0 ≤ δ < 1.

Run with Python 3 and SymPy 1.14.0:

    python3 docs/reports/michelson-trace-violation-two-dim/check.py

The checked output is in `output.txt`. All evidence uses exact symbolic algebra
or rational arithmetic. The program reproduces the squared-quotient limit in
Proposition 5.43 and evaluates the original principal-root definitions for a
noncommuting positive definite pair. It checks the determinant identities in the
proposed uniform upper-bound argument, all 2,800 cases of the displayed scalar
grid, 756 complex-vector cases, unitary transport, and the symbolic traces,
determinants and squared-quotient limit of a strictly positive sharpness family.

The finite grids are supporting tests. They do not prove the conjecture for all
matrices. The proposed settlement uses a general two-component inequality,
spectral decomposition and scalar Cauchy–Schwarz, followed by a sharpness limit.
No Lean proof is included.

Source: https://arxiv.org/abs/2610.34388v1
