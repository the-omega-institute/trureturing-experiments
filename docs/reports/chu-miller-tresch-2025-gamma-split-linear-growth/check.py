#!/usr/bin/env python3
"""Checks for H. V. Chu, S. J. Miller and G. Tresch, "Problems Regarding a Pair of Diophantine Equations",
arXiv:2512.12681v1, Section 3.4 (slower growth rate for given splits).

Gamma(a, b) = 0 iff (a/g) x + (b/g) y = (a/g - 1)(b/g - 1)/2 has a solution in nonnegative integers (g = gcd(a, b)), and
Gamma(a, b) = 1 iff 1 + (a/g) x + (b/g) y = (a/g - 1)(b/g - 1)/2 has one (Section 1; Theorem 1.1: exactly one does).
Fidelity: Theorem 1.1 (exactly one equation is solvable) for coprime pairs up to 80; Proposition 2.5 (Gamma(a, 2a) = 0,
Gamma(a, 2a - 1) = 1); Gamma(n, n + 1) = 0, 1, 0, 1, ... (quoted from [1, Theorem 1.6]); Theorem 2.1's inverse criterion.
Claim: the blocks Z = (3,7,17,19,23,43), O = (3,23,27,29,33,43) shifted by 40k have constant labels 0 and 1 (symbolic
certificates in k, and brute force for k <= 300); the scheduled sequence a_{5k+j} = 40k + B_k[j] with B_k = Z iff
floor((k+1)p) - floor(kp) = 1 is strictly increasing, satisfies -9 <= a_n - 8n <= 15, and 0 <= pN - C_N < 5 for N <= 4000
and several p. Exits 1 on failure.
"""
import math, sys
from fractions import Fraction as Fr
import sympy as sp

ok = True
def solvable(a, b, M):
    if M < 0: return False
    return any((M - a * x) % b == 0 for x in range(M // a + 1))
def gamma(a, b):
    g = math.gcd(a, b); a, b = a // g, b // g
    M = (a - 1) * (b - 1) // 2
    return 0 if solvable(a, b, M) else 1
def eq1(a, b):
    g = math.gcd(a, b); a, b = a // g, b // g
    return solvable(a, b, (a - 1) * (b - 1) // 2 - 1)
exactly_one = all((solvable(a, b, (a - 1) * (b - 1) // 2)) != eq1(a, b) for a in range(1, 81) for b in range(1, 81) if math.gcd(a, b) == 1)
print(f"fidelity Theorem 1.1: for coprime a, b <= 80 exactly one of the two equations is solvable: {exactly_one}")
p25 = all(gamma(a, 2 * a) == 0 for a in range(1, 301)) and all(gamma(a, 2 * a - 1) == 1 for a in range(2, 301))
print(f"fidelity Proposition 2.5 (a <= 300): Gamma(a, 2a) = 0 and Gamma(a, 2a - 1) = 1: {p25}")
alt = [gamma(n, n + 1) for n in range(1, 41)]
alt_ok = alt == [0, 1] * 20
print(f"fidelity Gamma(n, n + 1), n = 1..40, alternates 0, 1, 0, 1, ...: {alt_ok}")
def theta(a, b):
    g = math.gcd(a, b); return pow(a // g, -1, b // g)
t21 = True
for a in range(1, 61):
    for b in range(1, 61):
        if a % b == 0 or b % a == 0:
            t21 &= gamma(a, b) == 0; continue
        g = math.gcd(a, b)
        want0 = (theta(b, a) % 2 == 1) if (a // g) % 2 else (theta(a, b) % 2 == 1)
        t21 &= (gamma(a, b) == 0) == want0
print(f"fidelity Theorem 2.1 (a, b <= 60): {t21}")
ok &= exactly_one and p25 and alt_ok and t21
k = sp.symbols('k', integer=True, nonnegative=True)
rows = [(0,3,7,15*k+2,5*k),(0,7,17,6*k+2,14*k+2),(0,17,19,10*k+4,10*k+4),(0,19,23,15*k+8,5*k+2),(0,23,43,7*k+7,13*k+7),
        (1,3,23,13*k+7,7*k),(1,23,27,5*k+3,15*k+8),(1,27,29,10*k+7,10*k+6),(1,29,33,15*k+12,5*k+3),(1,33,43,6*k+6,14*k+11)]
cert = all(sp.expand(d + (40*k + c) * x + (40*k + e) * y - (40*k + c - 1) * (40*k + e - 1) / 2) == 0 and math.gcd(c, e - c) == 1 and 40 % (e - c) == 0 for d, c, e, x, y in rows)
print(f"certificates: the ten rows satisfy delta + a x + b y = (a - 1)(b - 1)/2 identically in k, with gcd(c, e - c) = 1 and (e - c) | 40: {cert}")
Z = (3, 7, 17, 19, 23, 43); O = (3, 23, 27, 29, 33, 43)
brute = all(gamma(40*kk + Z[j], 40*kk + Z[j+1]) == 0 and gamma(40*kk + O[j], 40*kk + O[j+1]) == 1 for kk in range(0, 301) for j in range(5))
print(f"blocks: Gamma = 0 on every edge of 40k + Z and Gamma = 1 on every edge of 40k + O for k <= 300 (brute force): {brute}")
ok &= cert and brute
def seq(p, n_terms):
    a = []
    K = n_terms // 5 + 2
    for kk in range(K):
        w = math.floor((kk + 1) * p) - math.floor(kk * p)
        B = Z if w == 1 else O
        a += [40 * kk + B[j] for j in range(5)]
    return a[:n_terms]
ps = [Fr(0), Fr(1), Fr(1, 2), Fr(1, 3), Fr(5, 7), math.sqrt(2) - 1, math.pi - 3, 0.123456]
NMAX = 4000
lab = {}
def G(a, b):
    if (a, b) not in lab: lab[(a, b)] = gamma(a, b)
    return lab[(a, b)]
for p in ps:
    a = seq(p, NMAX + 1)
    inc = all(a[i] < a[i + 1] for i in range(NMAX))
    growth = all(-9 <= a[n] - 8 * n <= 15 for n in range(NMAX + 1))
    C = 0; disc_ok = True; worst = 0.0
    for N in range(1, NMAX + 1):
        C += G(a[N - 1], a[N]) == 0
        d = float(p) * N - C
        worst = max(worst, abs(d))
        disc_ok &= -1e-9 <= d < 5
    print(f"p = {float(p):.6f}: strictly increasing {inc}; -9 <= a_n - 8n <= 15 {growth}; 0 <= pN - C_N < 5 for N <= {NMAX} {disc_ok}; C_N/N at N = {NMAX}: {C / NMAX:.6f}")
    ok &= inc and growth and disc_ok
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
