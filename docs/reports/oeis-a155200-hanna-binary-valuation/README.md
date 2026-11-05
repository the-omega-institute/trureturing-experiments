# Hanna's exact binary valuation conjecture (OEIS A155200)

Source: [OEIS A155200](https://oeis.org/A155200/internal), internal revision #39,
3 August 2026; valuation conjecture by Paul D. Hanna, 1 September 2009.
The ordinary generating function is
`A(x) = exp(sum_{k>=1} 2^(k^2) x^k/k)`.
The conjecture asserts `v_2([x^n] A(x)) = s_2(n)` for every `n >= 0`.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/oeis-a155200-hanna-binary-valuation/check.py
```

`output.txt` contains the output. The checker reproduces the entry's first
14 coefficients and independently compares its recurrence with direct truncated
rational exponential expansion through degree 24. It checks the valuation,
odd integral quotient, and normalized coefficient `n! a(n)/2^n == 1 (mod 4)`
through index 128, including the zero index and binary power boundaries.
All evidence uses exact integers and `fractions.Fraction`.
These are finite checks, not a proof of the universal conjecture.

Published integrality is separate from the valuation question: see Peter Bala's
comment in A155200, Jaidee–Moss–Ward, arXiv:1809.09199v2, and Sela Fried,
arXiv:2607.24832v1, Theorem 43. The falling-factorial recurrence is also an
instance of the standard weighted-permutation recurrence recorded in Ahaan
Kallat, arXiv:2610.15849, Proposition “Fundamental recurrence”.
