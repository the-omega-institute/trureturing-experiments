# OEIS A107099: the half-iterate of x + 4x^3 modulo 3

`check.py` computes A with A(A(x)) = x + 4x^3 through the integral series B = (A - x)/2, compares with the entry's b-file (`b107099.txt`, copied from oeis.org), and checks the settlement of the open-problem preregistration in `the-omega-institute/trureturing`: [x^n]A = 0 (mod 3) unless n = 3^k, and [x^(3^k)]A = binom(1/2, k) (mod 3).

Run: `python3 check.py 243 b107099.txt` (exit 0, `RESULT PASS`; output in `out.txt`).
