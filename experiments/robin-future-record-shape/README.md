# A future record need not have a monotone normalized pressure profile

For the original source `N = 2520`, the program certifies its strict
own-price inventory and the later inventory `C_t = 5040` at `t = 11`
and `t = 113/10`. With the complete pressure `P`, the tangent cap `c_A`
and `q_A(t) = P(t) - P(A)`, it bounds the derivative of `q_A/c_A`
strictly below zero at the first clock and strictly above zero at the
second. It also retains the positive `q_A - H(t, B(t))` margin.

Run with Python 3 and its standard library:

```sh
python3 experiments/robin-future-record-shape/certificate.py \
  --output experiments/robin-future-record-shape/result.json
```

The interval operations use outward integer fixed-point rounding. The
logarithm bounds reuse `log_enclosure` from the existing
`robin-density-transport-3e47/certificate.py`, including its explicit
positive-series remainder. Every occupied last layer and first absent
layer is checked; the inequality `r_(p,1) < 1/(p log p)` pays all primes
at least 13.

The global all-integer future-record property is a separate published
input: Caveney--Nicolas--Sondow, arXiv:1112.6010v2, Theorem 5(i)'s
proof, printed pp.14--15. Its equation (24), unconditional bound (2),
and exact set `R` give `G(m) < G(2520)` for every integer `m > 2520`.
The program neither reruns Robin's finite computation nor replaces
that argument by a finite future scan.

The result excludes a monotonicity claim with no lower restriction on
the source size. It says nothing against an additional `N > 5040` or
eventual-size hypothesis. No RH, stronger Robin bound, literature
priority, or Lean certification is claimed.
