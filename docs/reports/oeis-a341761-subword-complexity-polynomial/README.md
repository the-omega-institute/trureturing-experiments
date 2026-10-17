# OEIS A341761: the leading coefficient of the subword complexity polynomial

S(n, x) sums the number of distinct nonempty factors over all words of length n on x letters. S. Guo conjectures that its coefficient of x^n is n(n+1)/2.

Grouping words by the letters used gives S(n, x) = Σ_j C(x, j)·F(n, j). Only j = n reaches degree n, and F(n, n) = n!·n(n+1)/2.

`check.py` (Python 3) computes S(n, x) by brute force and interpolates. It then checks:
- the published rows;
- the leading coefficient;
- F(n, n).

Run `python3 check.py 6`. Exit 0 and `RESULT PASS` mean every check holds.
