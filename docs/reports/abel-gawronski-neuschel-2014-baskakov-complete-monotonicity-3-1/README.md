# Powered Baskakov complete monotonicity

Exact arithmetic checks for Conjecture 3.1 of Ulrich Abel, Wolfgang Gawronski
and Thorsten Neuschel, *Complete Monotonicity and Zeros of Sums of Squared
Baskakov Functions*, [arXiv:1411.7945v1](https://arxiv.org/abs/1411.7945v1).
The journal citation is Applied Mathematics and Computation 258 (2015), 130–137,
[DOI:10.1016/j.amc.2015.01.062](https://doi.org/10.1016/j.amc.2015.01.062).

Run `python3 check.py > check.out` with Python 3 and SymPy 1.14.0.
The final line must be `ALL_CHECKS_OK`. All comparisons use exact integers,
rationals or algebraic expressions; no floating point is used as evidence.

The fidelity check recovers the source's finite squared-Baskakov identity for
n = 1, 2, 3 from an independently generated rational generating function.
The settlement checks verify the geometric reduction at alpha = 1, the rational
Laplace identities for r = 2, 4, 6, 8, the algebraic inputs to a negative kernel
value, and a negative alternating derivative of order 75 at x = 2 for r = 8.
Two exact Taylor/derivative recurrences agree. Symbolic scaling and eight samples
check the infinite family n = c, with the wrong-sign derivative at x = 2/c.

The negative-kernel trigonometric inequality, passage to wrong signs on every
tail, and continuity transport to an open alpha interval require the analytic
arguments in the preregistration. The program does not certify those limits,
a value of epsilon, or a Lean proof. The lower-power identities support the
surviving alpha = 1 cases r = 2, 4, 6; no classification for r >= 10 is claimed.
