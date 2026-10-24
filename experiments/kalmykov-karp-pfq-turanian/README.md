# Kalmykov–Karp 2017 Conjecture 1: shifted Turanian coefficients of pFq

Numerical companion to the refutation of the coefficient clause of Conjecture 1 in
S. I. Kalmykov and D. B. Karp, *Log-concavity and Turán-type inequalities for the
generalized hypergeometric function*, Analysis Math. 43 (2017), arXiv:1611.06724.

Requires Python 3 with mpmath. Run from this directory:

    python3 coeff_scan.py
    python3 random_control.py 3 3000
    python3 turanian_scan.py 0.01 7 10
    python3 turanian_scan.py 0.5 6 11
    python3 turanian_scan.py 1 7 11

The kernel-checked result (delta_1 < 0 for p = 0, b = (1/2)^q, mu = 1,
alpha = beta = 1/2 and every q >= 6) lives in the main repository at
D5/S3/Analytic/HypergeometricShiftTuran/Result.lean.
