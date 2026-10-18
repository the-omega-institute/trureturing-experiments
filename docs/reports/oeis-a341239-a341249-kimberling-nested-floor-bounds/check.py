"""OEIS A341239 and A341249 (C. Kimberling): a(n) = floor(r*floor(s*n)) with s = sqrt(2) and r = 1 + sqrt(2),
resp. r = 2 + sqrt(2). Conjectures: 1 < r s n - a(n) < 3 (A341239) and 1 < r s n - a(n) < 2 (A341249).

Proof checked here: write s n = m + f with m = floor(s n), 0 < f < 1. Since sqrt(2) m = 2n - sqrt(2) f,
- r = 1 + sqrt(2): r m = m + 2n - sqrt(2) f, so a(n) = m + 2n - c with c = 1 if sqrt(2) f < 1, c = 2 otherwise,
  and r s n - a(n) = f + c, in (1, 1 + 1/sqrt(2)) or (2 + 1/sqrt(2), 3): the A341239 bound holds.
- r = 2 + sqrt(2): r m = 2m + 2n - sqrt(2) f, a(n) = 2m + 2n - c, and r s n - a(n) = 2f + c, in (1, 1 + sqrt(2)) or
  (2 + sqrt(2), 4). At n = 2: f = 2 sqrt(2) - 2, c = 2, value 4 sqrt(2) - 2 = 3.657...: the A341249 bound fails.
The script evaluates both with 80-digit arithmetic for n <= N, checks the data and the decompositions.
"""
import sys
from mpmath import mp, sqrt, floor, mpf
mp.dps = 80
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
s = sqrt(2)
D39 = [2, 4, 9, 12, 16, 19, 21, 26, 28, 33, 36, 38, 43, 45, 50, 53, 57, 60, 62, 67, 70, 74, 77, 79]
D49 = [3, 6, 13, 17, 23, 27, 30, 37, 40, 47, 51, 54, 61, 64, 71, 75, 81, 85, 88, 95, 99, 105, 109]
def run(r, D, lo, hi, kind):
    a = [int(floor(r * floor(s * n))) for n in range(1, N + 1)]
    okd = a[:len(D)] == D
    vals = [r * s * n - a[n - 1] for n in range(1, N + 1)]
    okb = all(lo < v < hi for v in vals)
    okdec = True
    for n in range(1, N + 1):
        m = floor(s * n); f = s * n - m; c = 1 if s * f < 1 else 2
        want = (f + c) if kind == 1 else (2 * f + c)
        if abs(vals[n - 1] - want) > mpf(10) ** -60: okdec = False
    return okd, okb, okdec, vals
o1 = run(1 + s, D39, 1, 3, 1)
o2 = run(2 + s, D49, 1, 4, 2)
fail2 = [n for n in range(1, 40) if not (o2[3][n - 1] < 2)]
print(f"A341239: data {o1[0]}; 1 < value < 3 for n<={N}: {o1[1]}; value = f + c: {o1[2]}")
print(f"A341249: data {o2[0]}; value in (1,4) for n<={N}: {o2[1]}; value = 2f + c: {o2[2]}; n < 40 with value >= 2: {fail2[:8]}; value at n = 2: {float(o2[3][1])}")
good = all(o1[:3]) and all(o2[:3]) and fail2[:1] == [2]
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
