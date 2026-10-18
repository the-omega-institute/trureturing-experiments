"""OEIS A380680 (P. D. Hanna): A(x) = Sum_{n=-oo..+oo} x^(n^2) (2 - x^n)^n/(1 - 2x^n)^n = 1 + 2 Sum_{n>=1} (...).
Conjecture: for n > 0, a(n) == 2 (mod 4) iff n is a square, and 4 | a(n) for nonsquare n.

Refutation checked here: a(1) = 4, so 1 is a square with a(1) == 0 (mod 4); also a(2) = 6 == 2 with 2 nonsquare.
Corrected statement: a(m) == 2 (mod 4) iff m = 2k^2, else a(m) == 0 (mod 4), because
a(m) = 2 c(m) with c = [x^m] Sum_{n>=1} x^(n^2) (2 - x^n)^n (1 - 2x^n)^(-n) and, modulo 2, each term is x^(2n^2).
The script expands the series exactly, compares with the published data, checks a(1), a(2), and the corrected
statement for m <= M.
"""
import sys
from sympy import symbols, series, expand
M = int(sys.argv[1]) if len(sys.argv) > 1 else 120
x = symbols('x')
s = 1
n = 1
while n * n <= M:
    s += 2 * x ** (n * n) * series((2 - x ** n) ** n / (1 - 2 * x ** n) ** n, x, 0, M + 1).removeO()
    n += 1
s = expand(s)
a = [int(s.coeff(x, m)) for m in range(M + 1)]
DATA = [1, 4, 6, 12, 32, 48, 120, 192, 450, 784, 1704, 3072, 6624, 12288, 25536, 49404]
ok_data = a[:len(DATA)] == DATA
from math import isqrt
def twice_sq(m): return m % 2 == 0 and isqrt(m // 2) ** 2 == m // 2
ok_ref = a[1] == 4 and a[2] == 6
ok_corr = all((a[m] % 4 == 2) == twice_sq(m) and (a[m] % 4 in (0, 2)) for m in range(1, M + 1))
print(f"data: {ok_data}; a(1) = {a[1]}, a(2) = {a[2]} (refutation): {ok_ref}; corrected: a(m) = 2 mod 4 iff m = 2k^2, else 0 mod 4, m<={M}: {ok_corr}")
good = ok_data and ok_ref and ok_corr
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
