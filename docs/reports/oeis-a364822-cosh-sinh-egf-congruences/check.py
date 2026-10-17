"""OEIS A364822 (M. Tebni): e.g.f. cosh(x)/(1 - 2 sinh(x)).

Checked statements:
- C1: for p prime, p > 2, a(p) = 2 (mod p).
- C2: for n = 2^m, a(n) = 1 (mod n); in fact v_2(a(n) - 1) = v_2(n) + 2 for every even n >= 2.

Proof ingredients checked here:
- a(n) = sum_{r>=1} (2^(r-1)/r) * (n+1)! [x^(n+1)] sinh(x)^r, from cosh*sinh^k = (sinh^(k+1))'/(k+1), and
  2^r sinh^r = sum_j C(r,j)(-1)^j e^((r-2j)x); this gives
  a(n) = sum_{r=1..n+1} (1/(2r)) sum_j C(r,j)(-1)^j (r-2j)^(n+1)   (compared with a series expansion).
- C1: from cosh(x) (2 sinh x)^m = sum_j C(m,j)(-1)^j (e^((m-2j+1)x) + e^((m-2j-1)x))/2,
  a(n) = sum_{m=0..n} sum_j C(m,j)(-1)^j ((m-2j+1)^n + (m-2j-1)^n)/2 (compared with the first form);
  by Fermat, a(p) = sum_m sum_j C(m,j)(-1)^j (m-2j) (mod p), where only m = 1 contributes, giving 2.
- C2: for even n, the r = 1 term is 1, the r = 3 term is 3^n - 1, and for odd r >= 5 the number of
  j < r/2 with C(r,j) odd and r - 2j = +-3 (mod 8) is even (checked for r <= R).
"""
import sys
from math import comb, factorial
from fractions import Fraction
from sympy import symbols, cosh, sinh, series, primerange
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 256
R = int(sys.argv[2]) if len(sys.argv) > 2 else 100001
def a(n):
    s = Fraction(0)
    for r in range(1, n + 2):
        s += Fraction(sum(comb(r, j) * (-1) ** j * (r - 2 * j) ** (n + 1) for j in range(r + 1)), 2 * r)
    assert s.denominator == 1
    return int(s)
x = symbols('x')
ser = series(cosh(x) / (1 - 2 * sinh(x)), x, 0, 16).removeO()
ok_series = all(ser.coeff(x, n) * factorial(n) == a(n) for n in range(16))
DATA = [1, 2, 9, 56, 465, 4832, 60249, 876416, 14570145, 272502272, 5662834089, 129446475776]
ok_data = [a(n) for n in range(len(DATA))] == DATA
def b(n):
    t = sum(comb(m, j) * (-1) ** j * ((m - 2 * j + 1) ** n + (m - 2 * j - 1) ** n) for m in range(n + 1) for j in range(m + 1))
    assert t % 2 == 0
    return t // 2
ok_second = all(b(n) == a(n) for n in range(40))
ok_fermat_sum = all(sum(comb(m, j) * (-1) ** j * (m - 2 * j) for m in range(N + 1) for j in range(m + 1)) == 2 for N in range(1, 60))
A = {}
def aa(n):
    if n not in A: A[n] = a(n)
    return A[n]
c1 = all(aa(p) % p == 2 for p in primerange(3, NMAX + 1))
def v2(z):
    c = 0
    while z % 2 == 0: z //= 2; c += 1
    return c
c2 = all(v2(aa(n) - 1) == v2(n) + 2 for n in range(2, NMAX + 1, 2))
c2pow = all(aa(2 ** m) % (2 ** m) == 1 % (2 ** m) for m in range(0, 9) if 2 ** m <= NMAX)
lemma = all(sum(1 for j in range(0, (r + 1) // 2) if (j & ~r) == 0 and (r - 2 * j) % 8 in (3, 5)) % 2 == 0 for r in range(5, R + 1, 2))
print(f"series n<16: {ok_series}; data: {ok_data}; second form n<40: {ok_second}; reduced sum = 2: {ok_fermat_sum}")
print(f"C1 odd p<={NMAX}: {c1}; C2 v2(a(n)-1) = v2(n)+2 for even n<={NMAX}: {c2}; a(2^m) = 1 mod 2^m: {c2pow}; parity lemma r<={R}: {lemma}")
good = ok_series and ok_data and ok_second and ok_fermat_sum and c1 and c2 and c2pow and lemma
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
