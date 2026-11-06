# Hanna A232692 positivity counterexample checks

Run `python3 check.py` with Python 3; no external packages or network are needed.
`output.txt` is its complete output. A successful run exits 0 and ends with
`ALL_CHECKS_OK`. All arithmetic uses integers or `fractions.Fraction`.

The derivative-polynomial recurrence and an independently implemented polynomial
ODE recurrence agree at every degree 0 through 430. The first negative EGF
coefficient for `(n,k)=(3,8)` is at degree 430. The full integer
`430! [x^430] G` is printed in `output.txt`. All 301 published table values agree.
A separate implementation of the literal integral/exponential fixed-point
equation agrees through degree 14 for seven parameter pairs. For `m=1..6`, exact
EGF multiplication verifies `G(x,3m,8m)^m = G(mx,3,8)` through degree 430 and finds
a negative coefficient in each family member. The all-`m` statement needs the
symbolic uniqueness and scaling argument; finite checks alone do not prove it.
These checks are experimental evidence, not Lean proofs.

Sources: Paul D. Hanna, [OEIS A232692](https://oeis.org/A232692/internal), revision
19, `%N`, `%C CONJECTURE`, `%e`, and its linked b-file; the threshold is
[OEIS A047399](https://oeis.org/A047399/internal), published offset 1. Related
source examples are [A232690](https://oeis.org/A232690) and
[A232691](https://oeis.org/A232691). `b232692.txt` is an unchanged source-table
snapshot from https://oeis.org/A232692/b232692.txt, credited there to Paul D.
Hanna. Its SHA-256 is
`fa11950a3994656977b98cc4ece6b4780565ad446ed198cd0e03e5568facec87`.
The table and quoted sequence data are credited to The Online Encyclopedia of
Integer Sequences and made available under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), as stated in the
[OEIS End-User License Agreement](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
