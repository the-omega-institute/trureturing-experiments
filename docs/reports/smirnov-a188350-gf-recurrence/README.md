# Smirnov's A188350 generating function and recurrence

Run `python3 check.py` using Python 3 and the standard library. The exact output is in `output.txt`; its final line is `ALL_CHECKS_OK`.

The source is [OEIS A188350](https://oeis.org/A188350/internal), revision 16, with the two FORMULA conjectures attributed to Mikhail Smirnov on October 15, 2026. R. H. Hardin contributed the native sequence in 2011. The [parent table A188358](https://oeis.org/A188358/internal) fixes the alphabet convention. `b188350.txt` is Hardin's [published table for n=1..200](https://oeis.org/A188350/b188350.txt), reproduced with attribution under the [OEIS End-User License Agreement](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement) (CC BY-SA 4.0).

A word has n+1 letters in {-1,0,1}, n edges, and free endpoints. The program independently accumulates the distribution of the adjacent-product sum and final letter. It checks the complete distributions by literal enumeration at n=0..8, matches all 200 source terms, and compares all coefficients with the collapsed two-state transfer through n=200. Generalized-binomial expansions with exact rational coefficients check the generating function, the normalized square root and the denominator-cleared binding equation through degree 200. The stated recurrence is checked for n=5..200.

These are finite exact checks. They support the proposed all-length transfer and constant-term argument but do not prove either conjecture for every n. No Lean proof or probe is included.
