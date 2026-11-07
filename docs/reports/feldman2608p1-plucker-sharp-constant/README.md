# Exact checks for Feldman P1 over complex scalars

The target is the unweighted coordinate Plücker norm in David Victor Feldman,
arXiv:2608.00983v1, §9 P1. The proposed optimal constant for n ≥ 1 and 2 < p < ∞ is
C(n,p) = n^(n/2 − n/p).

Run `python3 check.py` with Python 3 and SymPy. The deterministic output is
`check-output.txt`; successful completion ends in `ALL_CHECKS_OK`.

The program independently implements the Leibniz determinant and uses exact
Gaussian-integer arithmetic for the input checks. It reproduces the source's
Cauchy–Binet identity, Gram–Hadamard bound, crude factorial bound and Sylvester
determinant examples. It checks the proposed bound at p = 2, 4, 6, 8 on all 94
binary matrices with 1 ≤ n ≤ 2 and n ≤ M ≤ 3, and on 384 Gaussian-integer matrices
with 1 ≤ n ≤ 4 and n ≤ M ≤ n+3, with random seed 260800983.

For n = 1, …, 8 the Fourier matrices are evaluated in Q[x]/Φ_n(x), with
conjugation x ↦ x^(n−1). The program checks their Gram matrices, determinant
absolute squares, and equality at p = 3, 4, 6, 8. It squares the nonnegative
p-th-power equality to avoid radicals. These checks also recover C(1,p) = 1 and
the known bilinear value C(2,p)^p = 2^(p−2).

Every calculation uses exact arithmetic. No sampled case establishes the
universal inequality for real p or the passage to infinitely many coordinates.
That proposed settlement requires the general analytic argument in the
preregistration. The program contains no imported third-party proof code.
