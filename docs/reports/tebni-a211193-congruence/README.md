# Tebni's A211193 congruence: exact finite checks

Run `python3 check.py` with Python 3; only the standard library is used.
`check.out` is its output, ending in `ALL_CHECKS_OK`.

The program computes the coefficients of
`exp(exp((1+x)*log(1+x))-1)` in exact rational arithmetic through index 200,
checks an independent integer differential recurrence, the signed Stirling
formula for the inner exponential, and OEIS's Bell/Lehmer–Comtet formula.
All 201 normalized coefficients match the supplied OEIS table. It checks
`n | a(n)-1` and the inner coefficient divisibility for `2 <= n <= 200`.
The divided-power transfer identity is checked through index 40, including
negative block weights. These are finite checks, not an all-index proof or
Lean kernel verification.

Source: [OEIS A211193](https://oeis.org/A211193), revision #49,
with the conjecture attributed to Mélika Tebni; the generating function and
initial terms are those of that entry. `b211193.txt` is the unmodified
[0..200 table](https://oeis.org/A211193/b211193.txt), attributed there to
Alois P. Heinz. OEIS data are supplied under the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
The check prints the table's SHA-256 to bind the reference data.
The independent triangle recurrence is from
[OEIS A008296](https://oeis.org/A008296).
