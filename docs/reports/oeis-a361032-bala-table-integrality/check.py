"""OEIS A361032 (P. Bala): T(n,k) = F(n)*(4k)!/(k! (k+n+1)!^3), F(n) = (4n+4)!/(8 (n+1)!).
The entry says this choice of F(n) "appears to produce integer values".

For n, k <= M the script checks that every T(n,k) is an integer, that
v2(T(n,k)) = 3*(s2(n+k+1) - 1) (s2 = binary digit sum), the per-prime-power floor
inequality floor(4a)+floor(4b) >= 3*floor(a+b) used for odd p (on rationals a, b with
denominators up to 60), and the published antidiagonals.
"""
import sys
from math import factorial
from fractions import Fraction

M = int(sys.argv[1]) if len(sys.argv) > 1 else 45
DATA = [3, 315, 9, 46200, 280, 280, 7882875, 17325, 3675, 17325, 1466593128, 1513512]

def v2(x): return (x & -x).bit_length() - 1
ok_int = ok_v2 = True
T = {}
for n in range(M + 1):
    F = Fraction(factorial(4 * n + 4), 8 * factorial(n + 1))
    for k in range(M + 1):
        t = F * Fraction(factorial(4 * k), factorial(k) * factorial(k + n + 1) ** 3)
        if t.denominator != 1: ok_int = False; continue
        T[n, k] = t.numerator
        if v2(t.numerator) != 3 * (bin(n + k + 1).count("1") - 1): ok_v2 = False
anti2 = []  # ascending antidiagonals: T(s,0), T(s-1,1), ..., T(0,s)
for s in range(5):
    for k in range(s + 1): anti2.append(T[s - k, k])
ok_data = anti2[:len(DATA)] == DATA
from math import floor
ok_floor = all(floor(4 * Fraction(a, d)) + floor(4 * Fraction(b, d)) >= 3 * floor(Fraction(a + b, d))
               for d in range(1, 61) for a in range(0, 2 * d) for b in range(0, 2 * d))
print(f"n,k <= {M}: all integral: {ok_int}; v2 = 3(s2(n+k+1)-1): {ok_v2}; floor inequality: {ok_floor}; data: {ok_data}")
good = ok_int and ok_v2 and ok_floor and ok_data
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
