# OEIS A372527: residues mod 4 of the solution of A(A(x)^2) = x A(x)/(1 - 2A(x))

`check.py` solves the equation by its division-free integer recursion, compares with the entry's b-file (`b372527.txt`, copied from oeis.org), and checks the settlement of the open-problem preregistration in `the-omega-institute/trureturing`: for n >= 2, a(n) = 2 (mod 4) iff n = 2^k + 1, else 0 (mod 4); a(1) = 1. It also checks the proof steps over F_2 (x solves the equation mod 2; C = (A - x)/2 satisfies C(x^2) = xC + x^3 mod 2 and equals sum x^(2^k+1) mod 2).

Run: `python3 check.py 160 b372527.txt` (exit 0, `RESULT PASS`; output in `out.txt`).
