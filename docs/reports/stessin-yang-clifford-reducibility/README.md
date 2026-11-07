# Stessin–Yang Problem 7.12(a): Clifford-family check

`check.py` uses exact Gaussian-integer sparse matrices and SymPy exact algebra. It verifies the Pauli characteristic polynomial from the source's worked example, the Clifford involution and anticommutation relations, zero trace for every nonidentity monomial, and the quartic annihilator for every group monomial in dimensions 8 and 16 (64 and 256 ordered monomials). It also factors the dimension-8 symbolic characteristic polynomials for representative subset sizes 0 through 6, obtaining fourth powers of quadratics or squares of quartics.

Run from this directory with `python3 check.py`. The final output line is `ALL_CHECKS_OK`.
