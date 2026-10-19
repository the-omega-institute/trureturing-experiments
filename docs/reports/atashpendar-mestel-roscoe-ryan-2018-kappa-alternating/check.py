#!/usr/bin/env python3
"""Atashpendar-Mestel-Roscoe-Ryan (arXiv:1807.11609v1), Conjecture IV.4: the alternating word minimises
kappa^2(x) = sum_{1<=r,s<=m} C(r+s-2, r-1) C(2m-r-s, m-r) [x_r = x_s] over binary words of length m.

(1) Exhaustive for m <= 14: the minimum is attained only by the two alternating words (consistent with the conjecture).
(2) The periodic words (011001)^k: exact integer kappa^2 by the direct double sum and by the path-cost recurrence;
    reports the smallest k with kappa^2((011001)^k) < kappa^2(a_{6k}) = 6k C(12k-2, 6k-1), and certifies k = 58.
(3) The Bernstein identity kappa^2 = ((2n+1) C(2n,n) + E_n)/2, E_n/C(2n,n) = (2n+1) int_0^1 p_n^2, checked exactly for
    small words; and E_n / C(2n, n) for (011001)^k approaching 815/819 (Fraction values reported).
Exit 0 and RESULT PASS when the counterexample is certified.
"""
import itertools, sys
from math import comb
from fractions import Fraction as F

def kappa(x):
    m = len(x); n = m - 1
    return sum(comb(i + j, i) * comb(2 * n - i - j, n - i) for i in range(m) for j in range(m) if x[i] == x[j])

def defect(x):  # m C(2n,n) - kappa(x) by the path recurrence (integers only)
    m = len(x); p = [0] * m; d = [0] * m
    for i in range(m):
        for j in range(m):
            if i == 0 and j == 0: p[0] = 1; continue
            paths = p[j] + (p[j - 1] if j else 0); dd = d[j] + (d[j - 1] if j else 0)
            p[j] = paths; d[j] = dd + paths * (int((i + j) % 2 == 0) - int(x[i] == x[j]))
    return d[-1], p[-1]

ok = True
for m in range(2, 15):
    vals = [(kappa(x), x) for x in itertools.product('01', repeat=m)]
    mn = min(v for v, _ in vals); arg = [x for v, x in vals if v == mn]
    a = tuple('01'[i % 2] for i in range(m))
    ok &= mn == m * comb(2 * m - 2, m - 1) and len(arg) == 2 and a in arg
first = None
for k in range(1, 61):
    x = '011001' * k; m = 6 * k
    D, P = defect(x)
    ok &= P == comb(2 * m - 2, m - 1)
    if D > 0 and first is None: first = k
    if k in (1, 2, 5, 10, 20, 40, 57, 58, 60):
        print(f"k={k:2d} m={m:3d} sign(alt - x)={'+' if D > 0 else '-' if D < 0 else '0'}  E/C - 1 = {F(-2 * D, comb(2 * m - 2, m - 1))}"[:120])
x58 = '011001' * 58
D58, _ = defect(x58)
K58 = kappa(x58); A58 = 348 * comb(694, 347)
ok &= (A58 - K58 == D58) and D58 > 0
print(f"smallest k with (011001)^k below the alternating word: {first}")
print(f"m=348: kappa^2(a) - kappa^2((011001)^58) = {D58}")
print(f"digits {len(str(D58))}; direct double sum agrees: {A58 - K58 == D58}")
# (3) Bernstein identity on random small words
import random; random.seed(1)
for trial in range(20):
    m = random.randint(2, 9); x = [random.randint(0, 1) for _ in range(m)]; n = m - 1
    eps = [(-1) ** b for b in x]
    # int_0^1 p_n^2 with p_n = sum eps_i B_i^n, exact via Beta integrals
    I = sum(F(eps[i] * eps[j] * comb(n, i) * comb(n, j), comb(2 * n, i + j) * (2 * n + 1)) for i in range(m) for j in range(m))
    E = I * (2 * n + 1) * comb(2 * n, n)
    ok &= F(kappa(x)) == ((2 * n + 1) * comb(2 * n, n) + E) / 2
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
