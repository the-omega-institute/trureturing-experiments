# Hanna A386667: exact dyadic vanishing checks

Run `python3 docs/reports/a386667-dyadic-vanishing/check.py` with Python 3.
The program uses only the standard library and exact `fractions.Fraction`
arithmetic. Its output is in `check.out`.

Source: Paul D. Hanna, [OEIS A386667](https://oeis.org/A386667), revision 23,
with the A083563 recurrence and the reversion formula in source formula (2.a).
The program reproduces the first ten displayed coefficients of both series,
verifies the original functional equation and two-sided reversion through
degree 128, and evaluates the literal quotient coefficient for all 21 pairs
with 2 <= n <= 7 and 0 <= k <= n-2. Each coefficient is zero. The Lagrange
transfer is independently cross-checked against inverse-power coefficients.
The excluded boundary n=2, k=1 has coefficient 24.

These are finite exact sanity checks, not a proof of the unbounded conjecture.
The reciprocal is taken only for a series with nonzero constant coefficient.
No floating point arithmetic, external packages, or network access is used.
