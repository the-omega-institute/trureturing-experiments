# Principal logarithms for Ad-invariant Finsler norms

Exact arithmetic checks for the unnumbered conjecture following Proposition 4.3 in Gabriel Larotonda, *Short curves in homogeneous spaces of the unitary group for unitarily invariant norms*, arXiv:2610.31535v1, Section 4.1. Source: <https://arxiv.org/abs/2610.31535v1>.

Run with Python 3 and SymPy:

    python3 check.py > check.out

Angles are rational multiples of π; squared lengths are given in units of π². There is no floating point arithmetic. The output ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce Section 4's SU(5) example: the scalar principal logarithm has eigenangles 2π/5, while a traceless logarithm has spectral norm at least 8π/5, attained by a single winding shift. They also check the diagonal Frobenius instance of Theorem 2.1 on a specified finite grid of logarithm shifts. These grid comparisons are finite checks, not a proof of the theorem.

For n = 2 and c = 4, the norm F(iH)² = Tr(H²) + 4(Tr H)² and H = π diag(4/5, 1/2) give principal squared length 153π²/20. The alternative H′ = π diag(−6/5, 1/2) has the same exponential and squared length 73π²/20. The program checks the explicit unitary curve exp(itH′), its endpoints, constant left logarithmic derivative and exact length. It does not compute the distance's exact value: a strictly shorter competitor suffices.

The symbolic checks cover the quadratic-form decomposition, a nontrivial unitary conjugation, the rank-one winding-shift identity, and parameterized witnesses for every n ≥ 2 and c > −1/n with c ≠ 0. There are 56 additional exact parameter checks. The positivity and triangle inequality for the full norm follow mathematically from its positive definite inner product; finite checks alone do not establish universal admissibility. For negative c the norm is rescaled by 1/√(1 + nc), ensuring domination of the operator norm. For c = 0, principal-logarithm minimality is the established Frobenius result. In dimension one every admissible member is a scalar multiple of the ordinary absolute value; the nonzero-c classification is stated only for n ≥ 2.

These are preregistration checks. No Lean proof of the conjecture or its negation is included.
