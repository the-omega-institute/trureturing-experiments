# Exact two-distance entropy continuity checks

Source: Michael G. Jabbour and Nilanjana Datta, “Tightening continuity bounds for
entropies and bounds on quantum capacities”, arXiv:2310.17329v3,
DOI 10.1109/JSAIT.2024.3469929. The target is Conjecture 1, equation (108),
with the function in equations (18)–(19). The classical fidelity result is
Theorem 3 and its saturating distributions in equations (53)–(55).

Run `python3 check.py`. Python 3.9 or newer and the standard library suffice.
The expected output is stored in `output.txt`.

Probabilities and logarithm coefficients are exact fractions. Rational
logarithms are represented by coefficients of logarithms of primes.
Equality follows from prime factorization; signs are decided by comparing
integer products after clearing coefficient denominators. No floating point
or approximate logarithm is used.

The program checks Theorem 3 on every unordered pair of probability vectors
in dimensions 2 through 5, with a common denominator from 1 through 6. It
also checks the identity

    h(e) + e f_d(1, v/e) = f_d(e, v)

and the source's saturating distributions at every feasible pair e=a/n,
v=b/n with 1 <= b <= a <= n <= 24 and dimension 2 through 10. It checks
normalization, exact total variation and local distances, normalized
orthogonal Jordan parts, entropy attainment, integer and noninteger ratios,
the endpoint e=1, and the separate equal-state extension f_d(0,0)=0.
Two negative controls omit the remainder and replace the support ceiling
by a floor; both must be detected.

These finite checks establish the listed identities and classical cases.
They do not prove the conjecture for arbitrary noncommuting quantum states.
The proposed general argument uses Theorem 1 of Audenaert, Bergh, Datta,
Jabbour, Capel and Gondolf, arXiv:2408.15306v4, together with the original
classical theorem and the scaling identity. No Lean proof is included.
