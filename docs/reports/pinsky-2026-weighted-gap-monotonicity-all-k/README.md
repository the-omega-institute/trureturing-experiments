# Exact checks for Pinsky's weighted-gap question

Source: Ross G. Pinsky, *A new look at some aspects of one-dimensional random
sequential adsorption and its continuum limit*, arXiv:2608.04730v1, Theorem 2,
equation (1.3), and the question immediately after Table 1.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`output.txt` contains the output of that command. All evidence uses integers,
`fractions.Fraction`, or exact SymPy algebra; there is no floating point or
numerical quadrature.

The fidelity check reproduces the source's longest-gap identity at k = 2,
g_(2;1) = exp(-2), by exact symbolic integration. It also checks the polynomial
derivative identity behind the source's general longest-gap formula for
2 <= k <= 60. The source's erroneous Table 1 entry at g_(10;5) is excluded by
an independent rational enclosure, agreeing with Machina PM9/PM18.

For the inequality check, put P_k(s) = 2 sum_(1<=j<k) s^j/j and write
exp(P_k(s)) = sum_n c_n s^n. The nonnegative coefficients satisfy
n c_n = 2 sum_(1<=j<k, j<=n) c_(n-j). The code encloses these coefficients
by directed rounding on a common rational grid with denominator 10^60.
The total coefficient mass exp(P_k(1)) is enclosed by a Taylor polynomial
and a geometric bound on its positive remainder. After truncation at n = 2000,
the remaining mass is divided by the smallest omitted denominator
(n+l+1)(n+l+2). This bounds the beta integrals of the tail, then division
by exp(P_k(1)) gives rigorous rational bounds for g_(k;l).

All 1770 adjacent weighted-gap inequalities for 2 <= k <= 60 hold. Their
minimum slack occurs at (k,l) = (60,58) and lies in
[5352695183241/10^18, 5352695183242/10^18]. The proposed stronger bound
Delta_l >= g_(k;l)/(2k-1) also passes throughout the window; its residual
is at least 4592321596501/10^18.

A separate symbolic check verifies the polynomial divisions and factorial
moments used in a scalar bound for the proposed uniform analytic proof.
It computes 27728403640961/27912134250000, which is less than the nominated
414358946399/416821204800, itself strictly less than 1. The analytic
inequalities connecting these certificates to the universal claim are
proposed in the preregistration, rather than certified by this program.
The finite window is evidence for the proposed settlement, not a proof
for arbitrary k. No Lean proof is included.
