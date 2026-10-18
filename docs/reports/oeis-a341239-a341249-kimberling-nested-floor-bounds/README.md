# OEIS A341239 and A341249: Kimberling's bounds for floor(r·floor(√2 n))

With s = √2 and a(n) = ⌊r⌊sn⌋⌋, C. Kimberling conjectures two bounds:
- 1 < rsn − a(n) < 3 for r = 1 + √2 (A341239);
- 1 < rsn − a(n) < 2 for r = 2 + √2 (A341249).

Writing sn = m + f gives rsn − a(n) = f + c (A341239) and 2f + c (A341249), with c ∈ {1, 2} set by whether √2 f < 1. So:
- the A341239 bound holds;
- the A341249 values fill (1, 1+√2) ∪ (2+√2, 4), and the bound fails, first at n = 2 (value 4√2 − 2 ≈ 3.657).

`check.py` (Python 3 with mpmath, 80 digits) checks the data, the bounds and the decompositions for n ≤ N.

Run `python3 check.py 20000`. Exit 0 and `RESULT PASS` mean every check holds.
