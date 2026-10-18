# OEIS A379202: parity of the solution of 1/x = Sum_{n in Z} A^n (A^n + 2)^(n+1)

`check.py` folds the negative summands, checks that S(y) = Sum_{n in Z} y^n (y^n + 2)^(n+1) equals y^(-1) T with T an integer series congruent to 1 mod 2, computes A as the compositional inverse of B = 1/S (the entry's formula (5)), compares with the entry's b-file (`b379202.txt`, copied from oeis.org), checks equation (1) and that a(n) is even for 1 < n <= N, the settlement of the open-problem preregistration in `the-omega-institute/trureturing`.

Run: `python3 check.py 160 b379202.txt` with sympy 1.14.0 (exit 0, `RESULT PASS`; output in `out.txt`).
