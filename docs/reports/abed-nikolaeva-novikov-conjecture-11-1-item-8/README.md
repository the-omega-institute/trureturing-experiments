# Item 8: two prescribed contrasts

Exact checks for the final clause of Abed–Nikolaeva–Novikov, Conjecture 11.1,
arXiv:2610.34388v1, and the lower-bound construction of Proposition 5.22.

Run `python3 check.py` with Python 3 and SymPy. The output is in `check-output.txt`.
The program uses exact integers, fractions and symbolic algebra throughout.

It reproduces the published swap-block construction, its two-parameter formula,
and repeated even-dimensional blocks including both operator-contrast endpoints.
It checks the two scalar identities in the proposed orbit comparison, all small
integer spectral pairings specified in its output, and independently oriented
rational complex matrices. The latter checks use the exact positive square-root
formula for a two-by-two positive semidefinite matrix and compare squared
nonnegative quantities without rounding.

These checks support an expected proof. They do not certify the universal claim
or replace the subsequent Lean proof. No Lean proof is included.
