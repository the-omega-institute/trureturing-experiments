# Minimal-matroid magic-index exact checks

This report supports a preregistration for Masato Konoike, “On the magic
positivity of Ehrhart Polynomials of dilated polytopes”, arXiv:2504.21395v1,
Question 5.5. The target is m-index(B(T_{k,n})) = n-k for all n >= 2 and
1 <= k <= floor(n/2). Chen–Crigler–Ding–Paquin–Schennach,
arXiv:2610.32050v1, ask for the same full answer.

Run `python3 check.py`. It uses only Python 3's standard library, integers
and exact `fractions.Fraction` arithmetic. The stored output is `check.out`.
There is no floating point arithmetic or numerical root finding.

The source's Proposition 5.4 says D_{2,4} is not magic positive and proves
the thresholds in ranks 1 and 2. The program reproduces those results for
the printed ranges. It also independently counts lattice points satisfying
Ferroni's base-polytope inequalities for n <= 7 and dilation t <= 3.

The displayed D_{2,4} expansion in the source's rank-two proof contains
an extra X^3 term. The source enumerator gives magic coefficients
(1, -5/6, 1/6, 0), as verified by reconstruction. In particular D_{2,4}(-1)=0,
whereas that displayed expression would give -1. Its stated failure of
magic positivity and the proposed Question 5.5 settlement are unaffected.

The threshold is checked at every positive integer dilation through n-k+2
for n <= 28. For n <= 60 the program checks exact alternating evaluations
of the truncated binomial factor, its degree and leading coefficient,
D(0), D(1), D(-(n-k)), positivity of S((n-k)X), and the three rational
dilations n-k-1/2, n-k, n-k+1/2. Every magic-basis conversion is checked by
independent forward expansion.

These are finite checks. They do not prove the unbounded assertion, nor
certify a Lean proof. The proposed general argument uses alternating signs
and the intermediate value theorem to factor S into positive linear factors
with roots in (-k,-1), and the factor X+(n-k) for necessity.

Sources:
- https://arxiv.org/abs/2504.21395v1 (Question 5.5, Proposition 5.4, Section 5.2).
- https://arxiv.org/abs/2610.32050v1 (Question 1 and Conclusion).
- https://arxiv.org/abs/2003.02679v3 (Ferroni's enumerator and lattice inequalities).
