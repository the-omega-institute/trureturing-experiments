# Sitnik's problem 11 on the exponential remainder

Numerical companion to the refutation of the conjecture in problem 11 of
S. M. Sitnik, *Inequalities for the exponential remainder of the Taylor series*
(1993 preprint, arXiv:1608.08784): the forward differences in n of
g_n(x) = R_{n-1}(x)/R_n(x) are not nonnegative for every order k.

Requires Python 3 with mpmath. Run from this directory:

    python3 differences_scan.py

The kernel-checked witness (k = 3, n = 1, x = 1) lives in the main repository at
D5/S3/Analytic/ExpRemainderDifferences.lean.
