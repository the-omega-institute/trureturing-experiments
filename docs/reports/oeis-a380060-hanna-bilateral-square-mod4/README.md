# OEIS A380060: residues mod 4 of the bilateral series A = Sum_{n in Z} x^(n^2) (A - x^n)^n / (1 - x^n A)^n

`check.py` checks that the summands n and -n coincide, solves the equation by its division-free recursion, compares with the entry's b-file (`b380060.txt`, copied from oeis.org), and checks the settlement of the open-problem preregistration in `the-omega-institute/trureturing`: for n > 0, a(n) = 2 (mod 4) iff n is a square, else 0 (mod 4); with A = 1 + 2C, C = sum x^(m^2) (mod 2).

Run: `python3 check.py 400 b380060.txt` with sympy 1.14.0 (exit 0, `RESULT PASS`; output in `out.txt`).
