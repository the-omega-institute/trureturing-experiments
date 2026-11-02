# Bouafia's diversity-maximizing algorithm

Source: Philippe Bouafia, *Magnitude and diversity of trees*, arXiv:2605.03681v1,
Section 4.1, Algorithm 1 and the open problem after Theorem 4.
https://arxiv.org/abs/2605.03681v1

Run with Python 3 and SymPy installed:

```sh
python3 check.py
```

`check.out` is the resulting output. All certificate arithmetic uses integers,
`fractions.Fraction`, and exact symbolic algebra; no floating point is used.
The grid denominator is 10^40. Square roots are enclosed by integer square-root
inequalities. Exponentials are enclosed by 40-term Taylor sums and the remainder
bound in Mathlib's `Real.exp_bound'`, applied to distance/4 in [0,1], followed
by inversion and fourth powers. Interval Gaussian elimination certifies positive
pivots, invertibility and signs of both linear solves.

The fidelity check reproduces the source's symbolic two-point inverse and
magnitude formula. The main check certifies the nine exact terminating-decimal
points written in `check.py`. The intended algorithm retains indices
{1,3,4,5,6,8,9}, but the discarded seventh coordinate has (Zw)_7 < 1. An explicit
competitor with epsilon = 1/100000 has strictly smaller average similarity.
The literal printed stopping rule returns on the first solve and gives negative
probabilities at points 2 and 7. Thus neither reading is correct on this instance.

The same strict signs and failure hold uniformly when each coordinate is changed
by at most 1/10^12 and the configuration is scaled by any factor in
[1-1/10^10, 1+1/10^10]. This box check supplies a family and identifies the
mechanism: simultaneous deletion can discard a point that must re-enter.
These exact-arithmetic checks are independent computational evidence. No Lean
proof of the assertion or its negation is included.

Huntsman's arXiv:2509.19565v3, Appendix G.3, gives a related failure for the
scale-zero distance-matrix recursion. Its published theorem does not state the
fixed positive-scale exponential-kernel assertion checked here; novelty of a
broader failure mechanism is not claimed.
