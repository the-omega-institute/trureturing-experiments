# Step-m Eulerian parity positivity

This exact-arithmetic check supports a proposed settlement of Xin Chen and
Yudian Chen's parity-positivity conjecture. The primary sources are
[OEIS A392664](https://oeis.org/A392664),
[A392665](https://oeis.org/A392665), and
[A392666](https://oeis.org/A392666), which fix `m = -3, -4, -5` respectively.
The full integer range `m <= -3` is Conjecture 2 in their
[executive summary, second Zenodo version](https://doi.org/10.5281/zenodo.18452218)
(publication date 2026-02-01; PDF headed January 2026, CC BY 4.0).

Run `python3 check.py`. Python's standard library is sufficient. All arithmetic
used as evidence is exact integer arithmetic. The stored output is `output.txt`;
a successful run exits 0 and ends with `ALL_CHECKS_OK`.

Fidelity checks reproduce the 55 published terms through row 10 of each OEIS
entry and the source's nonnegativity checks for every even row through 400 at
those three parameters. Additional checks use `r = -m = 3,...,10` through row
403 and `r = 1000000` through row 30: 3,254 rows and 651,713 entries. They test
reflection, strict positivity on even rows, strict negativity of odd-row
interiors, the source row-sum product and even alternating sum, and both
second-column identities used by the proposed induction. A separate recursive
implementation cross-checks 702 entries through row 12.

For `b_n = T^(-r)(n,2)` and `n >= 2`, the checked identities are

```text
r*b_n = 2*(1 + (-1)^n*(r-1)^(n-1)) - r*(n-1)
b_(n+2) = (r-1)^2*b_n + r*(r-2)*n - r^2 + 2.
```

These are finite checks. They establish neither the unbounded conjecture nor
the absence of a prior proof. No Lean proof is part of this experiment.

The checker is an independent implementation of the displayed recurrence; it
does not copy the archived source toolkit. The embedded numerical source data
are attributed to the three OEIS entries above. Production: codex-cli research
flight `0d35a4278217dfadec726f42`, with prior disclosure `repo-prior-exposed`;
search candidate from flight `520ad2406ad056967eaa62f5`.
