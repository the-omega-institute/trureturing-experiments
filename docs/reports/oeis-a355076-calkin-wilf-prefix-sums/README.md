# OEIS A355076: integral prefix sums of the Calkin–Wilf enumeration

A355076(n) is the denominator of S(n) = Σ_{k≤n} fusc(k)/fusc(k+1). R. Sigrist conjectures a(n) = 1 for n = 2·4^k − 1 and n = 2·4^k − 2.

The ratios fusc(k)/fusc(k+1) list the Calkin–Wilf tree row by row, with these properties:
- row m has sum 3·2^(m−1) − 1/2 (children a/(a+b), (a+b)/b and the reciprocal symmetry of a row);
- hence S(2^(m+1) − 1) = (3·2^(m+1) − m − 4)/2, an integer iff m is even;
- the last term of the row is m + 1.

`check.py` (Python 3) checks, for 2^(m+1) ≤ N:
- the published data;
- the row sums and the reciprocal symmetry;
- the prefix formula and the last term;
- the conjecture.

Run `python3 check.py 65536`. Exit 0 and `RESULT PASS` mean every check holds.
