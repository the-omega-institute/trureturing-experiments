# Symmetric Gaussian feedback checks

Source: M.-J. So, J. Moran, Y. Lim, M.-S. Choi and H. Kwon,
“Local Gaussian bounds on the non-destructive discrimination of two-mode squeezed states”,
arXiv:2609.01373v1, Sections III.2–III.3 and Appendix B.

Run with Python 3.12 or later and SymPy 1.14:

    python3 docs/reports/smlck2609-glocc-symmetric/check.py

`check.out` is the program output. The program uses exact SymPy algebra and
outward-rounded rational interval arithmetic; no floating point is used as evidence.
It verifies the conditional Gaussian square completion, covariance conditioning,
centering, outcome-independent Schmidt parameter, Gaussian fidelity and the
source's symmetric squeezing formula for 120 exact cases, together with symbolic
identities for arbitrary parameters. Rational Taylor bounds enclose the stationary
point and score, reproduce the source's rounded 0.713804 and 0.626829, and check
weak-squeezing values at lambda=1/10 and 1/5 and the plateau at 1/2, 3/4 and 9/10.

These checks validate formulas and finite instances. They do not prove optimality
over feedback operations or the infinite-dimensional Schmidt overlap bound.
The positive settlement route is an analytic argument, not finite enumeration.
