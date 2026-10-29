# Hanna scaled reversion at q = 6 and q = 7

Run `python3 docs/reports/hanna-scaledrev-q67/check.py` from the repository root.
Python 3's standard library is sufficient. The program downloads both OEIS
b-tables and prints their SHA-256 digests; `output.txt` is its successful output.

Sources: Paul D. Hanna, [A393858](https://oeis.org/A393858), revision #10,
and [A393859](https://oeis.org/A393859), revision #8. Both entries attribute
the sequences to March 11, 2026. Their equation is
`A(x - x*A(q*x)/q) = x`, with `A(x) = sum(n>=1, a(n)*x^n)`.
The referenced tables contain indices 1 through 80. OEIS data attribution
and licensing: [OEIS licence agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The program solves this literal equation coefficient by coefficient. Put
`g[1]=1` and `g[n]=-q^(n-2)*a[n-1]` for n>=2. The coefficient at degree n
of `A(g)` is `a[n]+sum(1<=k<n, a[k]*[x^n]g^k)`. Setting it to zero
uses only integer additions and multiplications. The power coefficients
are computed by convolution. No conjectured parity or congruence is used
to produce the coefficients.

Checks:

- The first 80 integer coefficients match each published b-table exactly.
- A separate rational Horner composition verifies the literal divided
  source equation through degree 12.
- All first 500 positive-index coefficients are computed exactly in
  `Z/14Z` at q=6 and `Z/8Z` at q=7. The first computation checks both
  parity and congruence modulo 7; the second checks congruence modulo 8.
  These are images of the integer triangular recurrence, with no modular
  division. The first 80 residues also agree with the full integers.
- At q=6 both compositions with `X+X^2` are checked through degree 500
  over `F2`, using binomial coefficients and the square coefficient formula.

All checks pass and the final line is `ALL_CHECKS_OK`. These finite exact
checks support the preregistration; they prove none of the unbounded
conjectures. There is no Lean proof in this report.
