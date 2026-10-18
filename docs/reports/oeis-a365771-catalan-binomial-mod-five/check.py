"""OEIS A365771 (P. D. Hanna): a(n) = C(2n+1, n)/(2n+1) * C(3n-1, n). With A033042 = sums of distinct powers of 5
(base-5 digits in {0, 1}, including 0), the conjectures are:
(1) a(5m) = 4 (mod 5) and (2) a(5m+1) = 2 (mod 5) for m in A033042, m > 0;
(3) a(n) = 0 (mod 5) for n > 0 unless n or n - 1 equals 5m with m in A033042.

Proof checked here: for n >= 1, a(n) = 2 M(n) / (3 (n+1)) with M(n) = (3n)!/(n!)^3 = C(3n,n) C(2n,n).
- By Kummer (Legendre), v_5(M(n)) is the total carry of n + n + n in base 5. If n has a digit >= 2 then
  v_5(M(n)) > v_5(n+1) (t trailing 4s give carry >= 2t > t; otherwise v_5(n+1) = 0 < v_5(M(n))), so 5 | a(n).
- If all digits are 0 or 1, there is no carry, M(n) = 1 (mod 5) by Lucas, and a(n) = 2/(3(n+1)) (mod 5): 4 if
  n = 0 (mod 5) and 2 if n = 1 (mod 5).
The script checks the identity a(n) = 2M(n)/(3(n+1)), the valuation comparison, M(n) mod 5 and the conjectures.
"""
import sys
from math import comb
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
def a(n): return comb(2 * n + 1, n) // (2 * n + 1) * comb(3 * n - 1, n) if n > 0 else 1
def d01(m):
    while m:
        if m % 5 > 1: return False
        m //= 5
    return True
def v5(x):
    c = 0
    while x % 5 == 0: x //= 5; c += 1
    return c
DATA = [1, 2, 20, 280, 4620, 84084, 1633632, 33256080, 701149020, 15191562100]
ok_data = [a(n) for n in range(len(DATA))] == DATA
ok_id = all(3 * (n + 1) * a(n) == 2 * comb(3 * n, n) * comb(2 * n, n) for n in range(1, 3000))
ok_val = True; ok_lucas = True; c1 = c2 = c3 = True
for n in range(1, N + 1):
    M = comb(3 * n, n) * comb(2 * n, n)
    if d01(n):
        if M % 5 != 1: ok_lucas = False
    else:
        if not v5(M) > v5(n + 1): ok_val = False
    v = a(n) % 5
    if n % 5 == 0 and d01(n // 5) and v != 4: c1 = False
    if n % 5 == 1 and n > 1 and d01((n - 1) // 5) and v != 2: c2 = False
    excl = (n % 5 == 0 and d01(n // 5)) or ((n - 1) % 5 == 0 and d01((n - 1) // 5))
    if not excl and v != 0: c3 = False
print(f"data: {ok_data}; identity 3(n+1)a(n) = 2M(n), 1<=n<3000: {ok_id}; v5(M) > v5(n+1) when a digit >= 2: {ok_val}; M = 1 mod 5 for 0/1 digits: {ok_lucas}")
print(f"conjectures (1) {c1}, (2) {c2}, (3) {c3} for n <= {N}")
good = ok_data and ok_id and ok_val and ok_lucas and c1 and c2 and c3
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
