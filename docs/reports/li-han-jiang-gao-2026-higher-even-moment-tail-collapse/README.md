# Higher even moment tail collapse

This report checks the proposed all-order extension of Li–Han–Jiang–Gao, arXiv:2607.05226v2, Section 6(v).  For `m = 2, 3, 4`, the exact rational script evaluates the scalar extremizer with `t = 2, 3` and rational square parameter `a`, verifies its mass, mean, variance, and `2m`-th moment, and checks both one-sided and two-sided objectives.  It also enumerates every rational basic feasible solution of the finite moment LP on the listed support grid, with the moment inequality either slack or active.  The resulting finite LP optima equal the proposed value `p`; for `m = 2`, this is the source's tail formula `(K - 1) / ((t^2 - 1)^2 + K - 1)`.

All evidence uses Python `fractions.Fraction`; no floating-point arithmetic is used.  The global polynomial certificate `R_a(x^2) / R_a(t^2)` is recorded in the script and explains why the same value is an upper bound for unrestricted real laws.

Run:

```sh
python3 check.py > check.out
```

The final output line is `ALL_CHECKS_OK`.
