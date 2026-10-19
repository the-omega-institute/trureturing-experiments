#!/usr/bin/env python3
"""Roldan (arXiv:2604.09491v1), Conjecture 1.3: for distinct odd primes p, q and n = p^2 q^3, the divisor set
D* = {1, p^2, pq, q^2, p^2 q^2, p q^3} is the unique maximiser of the energy E(n, D) of the integral circulant graph
ICG(n, D) over all nonempty D of proper divisors of n.

(1) Direct: eigenvalue on the class gcd(j, n) = g is lambda_g = sum_{d in D} c_{n/d}(g) (Ramanujan sums, closed form
    mu(m/h) phi(m)/phi(m/h), h = gcd(g, m)); E = sum_g phi(n/g) |lambda_g|. All 2047 sets, for each (p, q) listed.
(2) The search seat's polynomial certificate: in u = p - 3, v = q - 3, every competitor's coefficientwise majorant U_X
    satisfies coeff(U_X) <= S - Delta, with S the optimum's exact coefficients (all 2048 Boolean matrices).
(3) A brute-force spectrum for (p, q) = (5, 3): eigenvalues of the 675 x 675 adjacency matrix of ICG(n, D*) agree with (1).
Exit 0 and RESULT PASS.
"""
import itertools, math, sys
import numpy as np

def factor(m):
    f = {}; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
def phi(m):
    r = m
    for p in factor(m): r = r // p * (p - 1)
    return r
def mu(m):
    f = factor(m)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)
def ram(m, j):
    h = math.gcd(j, m); return mu(m // h) * phi(m) // phi(m // h)

ok = True
def energy(p, q, D):
    n = p * p * q ** 3; divs = [p ** a * q ** b for a in range(3) for b in range(4)]
    return sum(phi(n // g) * abs(sum(ram(n // d, g) for d in D)) for g in divs)
for p, q in [(5, 3), (7, 3), (11, 3), (13, 3), (17, 3), (3, 5), (3, 7), (5, 7), (7, 5), (11, 5)]:
    n = p * p * q ** 3
    proper = [p ** a * q ** b for a in range(3) for b in range(4) if (a, b) != (2, 3)]
    Dstar = {1, p * p, p * q, q * q, p * p * q * q, p * q ** 3}
    best = None; E_star = energy(p, q, Dstar); second = -1
    for mask in range(1, 2 ** 11):
        D = {proper[i] for i in range(11) if mask >> i & 1}
        E = energy(p, q, D)
        if D != Dstar: second = max(second, E)
        if best is None or E > best[0]: best = (E, D)
    uniq = best[1] == Dstar and second < E_star
    ok &= uniq
    print(f"p={p} q={q} n={n} E(D*)={E_star} next={second} unique_max={uniq}")
# (2) polynomial certificate
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r
Arows = [[[0], [-1], [1]], [[-3, -1], [2, 1], [1]], [[6, 5, 1], [2, 1], [1]]]   # A(u)[a][c] as coefficient lists in u
Brows = [[[0], [0], [-1], [1]], [[0], [-3, -1], [2, 1], [1]], [[-9, -6, -1], [6, 5, 1], [2, 1], [1]],
         [[18, 21, 8, 1], [6, 5, 1], [2, 1], [1]]]
wA = [[6, 5, 1], [2, 1], [1]]; wB = [[18, 21, 8, 1], [6, 5, 1], [2, 1], [1]]
def U(X):
    tot = [[0] * 4 for _ in range(3)]
    for a in range(3):
        for b in range(4):
            L = [[0] * 4 for _ in range(3)]
            for c in range(3):
                for f in range(4):
                    if X[c][f]:
                        for i, x in enumerate(Arows[a][c]):
                            for j, y in enumerate(Brows[b][f]): L[i][j] += x * y
            C = [[abs(L[i][j]) for j in range(4)] for i in range(3)]
            for i in range(3):
                for j in range(4):
                    if C[i][j]:
                        for i2, x in enumerate(wA[a]):
                            for j2, y in enumerate(wB[b]):
                                assert i + i2 <= 2 and j + j2 <= 3
                                tot[i + i2][j + j2] += C[i][j] * x * y
    return tot
Xstar = [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]
S = [[1384, 1684, 678, 92], [1192, 1456, 588, 80], [266, 326, 132, 18]]
Delta = [[96, 136, 56, 8], [94, 128, 54, 8], [24, 32, 14, 2]]
ok &= (U(Xstar) == S)
cells = [(c, f) for c in range(3) for f in range(4) if (c, f) != (2, 3)]
cert = True
for mask in range(2 ** 11):
    X = [[0] * 4 for _ in range(3)]
    for i, (c, f) in enumerate(cells):
        if mask >> i & 1: X[c][f] = 1
    if X == Xstar: continue
    Ux = U(X)
    if any(Ux[i][j] > S[i][j] - Delta[i][j] for i in range(3) for j in range(4)): cert = False; print("cert fail", X); break
ok &= cert
# exact value check of the polynomial model against (1) at (p, q) = (5, 3): E(D*) = S(u=2, v=0)
ok &= sum(S[i][j] * 2 ** i * 0 ** j for i in range(3) for j in range(4)) == energy(5, 3, {1, 25, 15, 9, 225, 135})
print(f"certificate over 2047 competitors: {cert}; S equals U(X*): {U(Xstar) == S}")
# (3) brute-force spectrum for (5, 3)
p, q = 5, 3; n = p * p * q ** 3; Dstar = {1, 25, 15, 9, 225, 135}
Adj = np.zeros((n, n))
for a in range(n):
    for b in range(n):
        if a != b and math.gcd((a - b) % n, n) in Dstar: Adj[a, b] = 1
ev = np.linalg.eigvalsh(Adj)
ok &= abs(np.abs(ev).sum() - energy(5, 3, Dstar)) < 1e-6
print(f"brute-force energy (5,3): {np.abs(ev).sum():.6f} vs formula {energy(5, 3, Dstar)}")
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
