"""OEIS A367191 (J.-F. Alcover): a(n) is the largest positive k such that k^2 + k*n + n^2 is a perfect
square, or 0 if no such k exists. Conjectured formulas by n mod 4:
  n = 0 (mod 4):      (n/4 - 1)(3n/4 + 1),
  n = 1, 3 (mod 4):   4((n+3)/4 - 1)(3(n+3)/4 - 2) = (n-1)(3n+1)/4,
  n = 2 (mod 4):      2((n+2)/4 - 1)(3(n+2)/4 - 1) = (n-2)(3n+2)/8.

Proof checked here: with u = 2m - (2k+n) and v = 2m + (2k+n), the equation k^2 + kn + n^2 = m^2 becomes
u v = 3n^2, and k = (v - u - 2n)/4, m = (u + v)/4. k decreases as u grows, so the largest k comes from the
least admissible u: u = 1 (n odd), u = 2 (n = 2 mod 4), u = 4 (n = 0 mod 4).

The script (1) finds a(n) by brute force over 1 <= k <= (3n^2 - 2n - 1)/4 for 1 <= n <= N1 and compares
with the published data and the formulas; (2) for 1 <= n <= N2 checks the formula's witness k with the
witness m = (u + v)/4, and that no u smaller than the stated one is admissible (u v = 3n^2, u = v mod 2,
u + v = 0 mod 4, v - u - 2n = 0 mod 4, k > 0).
"""
import sys
from math import isqrt
from fractions import Fraction as Fr
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 150
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
DATA = [0,0,0,5,0,16,10,33,7,56,32,85,20,120,66,161,39,208,112,261,64,320]
def formula(n):
    if n % 4 == 0: return Fr(n, 4) - 1, Fr(3 * n, 4) + 1, 1
    if n % 4 == 2: return 2 * (Fr(n + 2, 4) - 1), Fr(3 * (n + 2), 4) - 1, 1
    return 4 * (Fr(n + 3, 4) - 1), Fr(3 * (n + 3), 4) - 2, 1
def fval(n):
    x, y, _ = formula(n); return x * y
def brute(n):
    best = 0
    for k in range(1, (3 * n * n) // 4 + 2):
        s = k * k + k * n + n * n; r = isqrt(s)
        if r * r == s: best = k
    return best
b = [brute(n) for n in range(0, N1 + 1)] if True else []
b[0] = 0
ok_data = b[:len(DATA)] == DATA
ok_brute = all(b[n] == fval(n) for n in range(1, N1 + 1))
print(f"brute force n<={N1}: matches data {ok_data}; matches formulas {ok_brute}")
def admissible(n, u):
    if (3 * n * n) % u: return False
    v = 3 * n * n // u
    return u < v and (u - v) % 2 == 0 and (u + v) % 4 == 0 and (v - u - 2 * n) % 4 == 0 and v - u - 2 * n > 0
ok_w = True
for n in range(1, N2 + 1):
    u0 = 1 if n % 2 else (2 if n % 4 == 2 else 4)
    k = fval(n)
    if k > 0:
        v = 3 * n * n // u0; m = (u0 + v) // 4
        if not (admissible(n, u0) and k == (v - u0 - 2 * n) // 4 and k * k + k * n + n * n == m * m): ok_w = False
    if any(admissible(n, u) for u in range(1, u0)): ok_w = False
print(f"witness and least admissible u for n<={N2}: {ok_w}")
good = ok_data and ok_brute and ok_w
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
