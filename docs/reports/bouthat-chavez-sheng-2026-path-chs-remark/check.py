"""Bouthat-Chavez-Sheng, arXiv:2606.17356v1, Remark 4.3: the CHS norm of the path P_n.
||P_n||_{X,d}^d = h_d(lambda(P_n)) = [t^d] 1/det(I - t A_n), computed exactly from the recurrence
det(I - tA_n) = det(I - tA_{n-1}) - t^2 det(I - tA_{n-2}).

Checks, for 1 <= n <= N and even 2 <= d <= D:
  (A1) the trigonometric expression (mpmath, 80 digits);
  (A2) the bilateral binomial sum (exact);
  (B)  the stated truncation of the bilateral sum to j in [-floor((d-2)/(2(n+1))), 1 + floor(d/(2(n+1)))];
  (C)  the restricted closed form (n^2+n-d)/((n+d)(n+d+1)) * binom(n+d+1, d/2) under n+1 >= d/2.
Usage: python3 check.py N D
"""
import sys
from fractions import Fraction
from math import comb
import mpmath as mp
mp.mp.dps = 80

def C(a, b): return comb(a, b) if 0 <= b <= a else 0

def h_path(n, d):
    D = [[1], [1]]
    for m in range(2, n + 1):
        a, b = D[-1], D[-2]
        c = a[:] + [0] * (len(b) + 2 - len(a))
        for i, x in enumerate(b): c[i + 2] -= x
        D.append(c)
    P = D[n]
    inv = [Fraction(0)] * (d + 1); inv[0] = Fraction(1, P[0])
    for k in range(1, d + 1):
        inv[k] = -sum(Fraction(P[i]) * inv[k - i] for i in range(1, min(k, len(P) - 1) + 1)) / P[0]
    return inv[d]

def bilateral(n, d, lo=None, hi=None):
    if lo is None: lo, hi = -d - 2, d + 2
    return sum(C(n + d - 1, (n + 1) * j + d // 2) - C(n + d - 1, (n + 1) * j + d // 2 - 1) for j in range(lo, hi + 1))

def trig(n, d):
    s = mp.mpf(0)
    for k in range(1, n + 1):
        x = mp.pi * k / (n + 1)
        s += (-1) ** (k + 1) * mp.cos(x) ** (n - 1 + d) * mp.sin(x) ** 2
    return mp.mpf(2) ** (n + d) / (n + 1) * s

def main():
    N, Dm = int(sys.argv[1]), int(sys.argv[2])
    a1 = a2 = 0; bad_a1, bad_a2, bad_b, bad_c, bound_cases = [], [], [], [], 0
    for n in range(1, N + 1):
        for d in range(2, Dm + 1, 2):
            h = h_path(n, d)
            if abs(trig(n, d) - h) > mp.mpf(10) ** -50 * max(1, abs(h)): bad_a1.append((n, d))
            if bilateral(n, d) != h: bad_a2.append((n, d))
            if bilateral(n, d, -((d - 2) // (2 * (n + 1))), 1 + d // (2 * (n + 1))) != h: bad_b.append((n, d))
            if n + 1 >= d // 2:
                r = Fraction(n * n + n - d, (n + d) * (n + d + 1)) * comb(n + d + 1, d // 2)
                if n + 1 == d // 2: bound_cases += 1
                if r != h: bad_c.append((n, d, h - r))
    print(f"(A1) trigonometric expression = h_d for n <= {N}, even d <= {Dm}: {bad_a1 == []}")
    print(f"(A2) bilateral sum = h_d: {bad_a2 == []}")
    print(f"(B) stated truncation fails in {len(bad_b)} cases; first {bad_b[:4]}")
    print(f"(C) restricted closed form fails in {len(bad_c)} cases, all with n+1 = d/2: {all(2 * (n + 1) == d for n, d, _ in bad_c)}; boundary cases checked {bound_cases}; differences h - formula: {sorted(set(x for _, _, x in bad_c))}; first {[(n, d) for n, d, _ in bad_c[:4]]}")
    print(f"witness n = 2, d = 6: h = {h_path(2, 6)}, formula = {Fraction(4 + 2 - 6, 8 * 9) * comb(9, 3)}")
    ok = bad_a1 == [] and bad_a2 == [] and bad_c != [] and all(2 * (n + 1) == d for n, d, _ in bad_c) and h_path(2, 6) == 1
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
