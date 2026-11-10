# Lewis–Souza proportional divisibility windows

For fixed rational `alpha` in `(0,1)`, let the window consist of the positive
integers `a` with `alpha*n < a <= n`, ordered by divisibility. The proposed
general argument for convergence of its order dimension uses a uniform finite
bound and dilation into every sufficiently large window.

Run with Python 3.12 or later, using only the standard library:

```sh
python3 docs/reports/lewis_souza_divisibility_window_limit/check.py
```

The stored output is `output.txt`. Every arithmetic comparison uses integers or
`fractions.Fraction`; no floating point, numerical optimization or random input
is used. Successful execution ends with `ALL_CHECKS_OK`.

The source fidelity checks reproduce the antichain result of Section 5 of
David Lewis and Victor Souza, *The order dimension of divisibility*,
arXiv:2001.08549v2, and the introduction's dimension formula for divisor sets
of `2, 4, 6, 12, 30`. The antichain dimension check requires at least two
elements. The dimension routine uses the literal minimum number of factors
in a real product, allowing the zero-fold product for a singleton. Assigning
singleton dimension one instead affects only finitely many initial windows
for each fixed alpha and does not affect the proposed limit assertion.

For every reduced rational `alpha=p/q` with `2 <= q <= 12` and `0 < p < q`,
the checks verify the small-prime coordinate embedding for `1 <= n <= 60`.
They also check the proposed dilation threshold for `1 <= n0 <= 30`, every
residue `0 <= s < n0`, and each of `q=Q, Q+1, 2Q+1`. Every source point is
checked for strict lower-endpoint membership and upper-endpoint membership;
every ordered pair is checked for preservation and reflection of divisibility.
Here `delta=floor(alpha*n0)+1-alpha*n0`, `Q=floor(alpha*n0/delta)+1`, and
the threshold is `Q*n0`.

Exact dimensions are checked for the declared smaller rational grid
`2 <= q <= 6`, `1 <= n <= 10`, whenever the window has at most eight
elements. All linear extensions are enumerated, and minimum realizer size is
an exact finite set-cover problem on ordered incomparable pairs. Chains and
antichains use explicit coordinate certificates and the matching elementary
lower bound. The output includes the dimensions and comparisons under the
guaranteed dilation threshold.

These are finite checks of the proposed route, not a proof of convergence
for all real alpha. No Lean proof is included.

Sources:

- David Lewis and Victor Souza, *The order dimension of divisibility*,
  JCTA 179 (2021), 105391, DOI 10.1016/j.jcta.2020.105391;
  <https://arxiv.org/abs/2001.08549>.
- Milan Haiman, *The Dimension of Divisibility Orders and Multiset Posets*,
  arXiv:2201.12952v5, Proposition 2.2 (the credited small-prime method).
- Victor Seixas Souza, *Combinatorial aspects of synchronisation,
  percolation, and arithmetic*, Cambridge PhD thesis, Question 6.10.3.
