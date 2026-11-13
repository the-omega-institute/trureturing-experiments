# ASM imprecise-copula permutation checks

`check.py` independently checks the proposed rounded-independence obstruction to
the unnumbered permutation-copula question in Section 8 of T. Košir, E. Perrone
and N. Stopar, *Dense alternating sign matrices and imprecise copulas*,
[arXiv:2607.06178v2](https://arxiv.org/abs/2607.06178v2).

Run with Python 3.12 or newer, using only the standard library:

```sh
python3 check.py
```

All evidence uses integers or `fractions.Fraction`. `check.out` is the complete
output; its final line is `ALL_CHECKS_OK`.

- Fidelity: reproduce the source's Examples 2.1 and 2.2, including scaling,
  forward and inverse corner sums, the ASM, and boundary-rectangle inequalities.
- For every order 1 through 40, check both ASMs, all four imprecise-copula
  inequalities on every rectangle (including degenerate ones), and the feasible
  real copula C(r,s)=rs/n, using multiplication by the positive denominator n.
- Exhaust all permutations at orders 1 through 8, reporting how many lie in
  the rounded bracket. These orders have witnesses.
- Check the dyadic Haar basis and its pairwise orthogonality at orders
  8, 16, 32 and 64. Verify empty-block counting, boundary cancellation, the
  coefficient identity and the integer sum-of-squares identity on all 40,320
  permutations of order 8 and 18 deterministic permutations at each larger
  order.
- Check the exact final contradiction arithmetic at exponents 2051, 2052 and
  4096 without materializing the grids.

The candidate family is P(r,s)=floor(rs/n), Q(r,s)=ceil(rs/n), n=2^m,
m≥2051. The general proposed argument uses at least n(m−2) empty dyadic
blocks, producing

```text
2048 * sum E² - n³ * t = 2 * sum (32 * E + n * H)² ≥ 0,
t ≥ n * (m - 2),
sum E² ≤ n⁴ under the permutation-bracket hypothesis.
```

These imply m≤2050. The program checks the identities and finite instances;
it does not prove the universal argument. It does not enumerate the
permutations of order 2^2051, certify the universal family, or provide a Lean
proof. No smallest counterexample order is asserted. The mathematical
argument and source fidelity remain subject to subsequent formal verification.
