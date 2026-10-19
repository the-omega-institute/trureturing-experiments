#!/usr/bin/env python3
"""OEIS A271960 (P. D. Hanna) conjecture: [x^(2^n - 2^k)] (x/A(x))^(2^n) = 0 for k = 0..n-2, n >= 2,
where A(x) = x + 2x^2 + ... satisfies A(x)^2 = A((x + 2 A(x)^2)^2).

1. Exact (integers): B = y - 2y^2 + (odd terms) from H(y)^2 = B(y^2), H = B + 2y^2 odd (the proof's Lemma 1
   structure), to order M = 520; A = series reversion of B by Lagrange inversion modulo two primes.
2. The defining equation A^2 = A((x + 2A^2)^2) is checked to order M modulo both primes, so the computed A is
   the sequence's g.f. (independently of the structure used in step 1); the OEIS EXAMPLE terms are checked exactly.
3. The conjecture: exactly (integer arithmetic) for 2^n <= 128, and modulo both primes for 2^n <= 512.
   The boundary k = n - 1 is reported (nonzero at n = 2, as the entry's range excludes it).
Exit 0 and RESULT PASS when all hold.
"""
import sys
import numpy as np
M = 520
PRIMES = [67108859, 67108837]

# 1. exact B
h = [0] * (M + 1); h[1] = 1
B = [0] * (M + 1); B[1] = 1; B[2] = -2
for m in range(2, M // 2 + 2):
    # coefficient of y^(2m) in H^2 equals B[m]
    if 2 * m - 1 > M: break
    s = sum(h[i] * h[2 * m - i] for i in range(3, 2 * m - 2, 2))
    bm = B[m] if m <= M else 0
    val = bm - s; assert val % 2 == 0
    h[2 * m - 1] = val // 2
    if 2 * m - 1 <= M: B[2 * m - 1] = h[2 * m - 1]
# B even coefficients beyond 2 are zero; B odd = h

def mulp(a, b, p, m):
    c = np.convolve(np.array(a[:m + 1], dtype=np.int64) % p, np.array(b[:m + 1], dtype=np.int64) % p)
    return list((c[:m + 1] % p).astype(np.int64))

def inv_unit(u, p, m):  # inverse of series with u[0] = 1 mod p
    r = [0] * (m + 1); r[0] = 1
    ua = np.array(u[:m + 1], dtype=object)
    for t in range(1, m + 1):
        r[t] = (-sum(int(u[i]) * r[t - i] for i in range(1, t + 1))) % p
    return r

ok = True
Amod = {}
for p in PRIMES:
    Bp = [x % p for x in B]
    yB = inv_unit(Bp[1:] + [0], p, M)            # y / B(y)
    A = [0] * (M + 1); P = [1] + [0] * M
    for n in range(1, M + 1):
        P = mulp(P, yB, p, M)                     # (y/B)^n
        A[n] = P[n - 1] * pow(n, -1, p) % p        # [x^n]A = (1/n)[y^(n-1)](y/B)^n
    Amod[p] = A
    # 2. defining equation
    A2 = mulp(A, A, p, M)
    inner = [(A2[t] * 2 + (1 if t == 1 else 0)) % p for t in range(M + 1)]
    g = mulp(inner, inner, p, M)
    comp = [0] * (M + 1); pw = [1] + [0] * M
    for i in range(1, M + 1):
        pw = mulp(pw, g, p, M)
        if A[i]:
            for t in range(M + 1): comp[t] = (comp[t] + A[i] * pw[t]) % p
    if comp != A2: ok = False; print("FAIL equation mod", p)
    # 3. modular conjecture
    xA = inv_unit(A[1:] + [0], p, M)
    for n in range(2, 10):
        N = 2 ** n
        Q = [1] + [0] * N
        base = xA[:N + 1]; e = N
        while e:
            if e & 1: Q = mulp(Q, base, p, N)
            base = mulp(base, base, p, N); e >>= 1
        for k in range(0, n - 1):
            if Q[N - 2 ** k] % p: ok = False; print("FAIL mod", p, n, k)
# exact check for small n
terms = [0, 1, 2, 9, 50, 312, 2086, 14613, 105864, 786627, 5962110, 45914544, 358247214]
for p in PRIMES:
    if [x % p for x in terms] != Amod[p][:13]: ok = False; print("FAIL terms mod", p)
def mul_exact(a, b, m):
    c = [0] * (m + 1)
    for i, x in enumerate(a[:m + 1]):
        if x:
            for j in range(0, m + 1 - i):
                if b[j]: c[i + j] += x * b[j]
    return c
# exact A to order 129 by Lagrange over the rationals (A has integer coefficients)
from fractions import Fraction
m = 129
u = B[1:m + 2]
yBe = [Fraction(0)] * (m + 1); yBe[0] = Fraction(1)
for t in range(1, m + 1):
    yBe[t] = -sum(u[i] * yBe[t - i] for i in range(1, t + 1))
Ae = [0] * (m + 1); P = [Fraction(1)] + [Fraction(0)] * m
for n in range(1, m + 1):
    P = mul_exact(P, yBe, m); v = P[n - 1] / n; assert v.denominator == 1; Ae[n] = int(v)
if Ae[:13] != terms: ok = False; print("FAIL exact terms")
ue = Ae[1:] + [0]
xAe = [0] * (m + 1); xAe[0] = 1
for t in range(1, m + 1):
    xAe[t] = -sum(ue[i] * xAe[t - i] for i in range(1, t + 1))
exact_checked = 0
for n in range(2, 8):
    N = 2 ** n; Q = [1] + [0] * N
    for _ in range(N): Q = mul_exact(Q, xAe, N)
    for k in range(0, n - 1):
        if Q[N - 2 ** k] != 0: ok = False; print("FAIL exact", n, k, Q[N - 2 ** k])
        exact_checked += 1
    print(f"n={n}: boundary k=n-1 coefficient = {Q[N - 2 ** (n - 1)]}")
print(f"order={M} exact_cases={exact_checked} (2^n<=128) modular_cases=36 (2^n<=512, two primes)")
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
