# AR(1) persistence and log-convexity

Exact checks for Titouan Donnart and Thomas Simon, *Persistence probabilities of
autoregressive chains with continuous innovations*, arXiv:2604.05670v1,
Remark 3(a).

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/autoregressive-persistence-logconvexity/check.py
```

The stored output is `check-output.txt`. All evidence uses integers and
`fractions.Fraction`; there is no floating point, random sampling, or numerical
quadrature.

The fidelity calculation directly convolves the killed symmetric bi-exponential
kernel at drift 1. Exact polynomial integration reproduces the source's
Section 4.2 formulas: p_n = binom(2n,n)/4^n and q_n = C_n/2^(2n+1), through
q_12. This is independent of the candidate's finite-state recurrence.

For the candidate, exact interval propagation classifies every innovation sign
word of lengths 0 through 12 on four pairs of bands, with drift 1/2, 1/3, 1/4,
and 1/5. These classifications agree with the independent recognizer “first sign
positive, no consecutive negatives”. Weighted sums for 27 rational positive-band
masses agree with the second-order recurrence. Its adjacent minors are checked
through index 64, its two-state killed-matrix representation through time 16,
and the source's supermultiplicativity inequality for indices 0 through 32.

The exact atomless witness can use uniform component laws on [12/5,13/5] and
[-11/10,-9/10], with positive-component mass 2/7. The one-negative lower margin
is 1/10 and the two-negative upper margin is -1/20. It has mean zero and
p_1=p_2=2/7, p_3=48/343, so p_1 p_3-p_2^2=-100/2401.

The script also checks the rational arithmetic in the Gaussian-contamination
estimate and the band-margin formulas for a parameter grid. It does not integrate
a Gaussian or verify smoothness. Finite checks do not establish the all-time
family theorem, the coupling bound, or the existence of smooth densities; these
require mathematical proofs. No Lean proof of the assertion or its negation is
provided here.
