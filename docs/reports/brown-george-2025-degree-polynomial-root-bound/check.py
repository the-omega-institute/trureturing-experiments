"""Check of Conjecture 4.1 of Brown and George, "On the roots of degree polynomials" (AJC 92(2), 2025; arXiv:2505.04882v1).

D(G; x) = sum_v x^deg(v).  For every graph on n <= 6 vertices (all edge sets) and 20000 random graphs on 7 vertices the
script computes the nonzero roots of D(G; x) and the ratio |z| / (n-1)^(|arg z|/pi); the conjecture says the ratio
is at most 1.  It also checks that the stars K_{1,n-1} attain the bound at z = (n-1)^(1/(n-2)) exp(i pi/(n-2)), and the
source's Theorem 3.2 (real roots: |z| <= n - 1) on the same graphs.
"""
import itertools
import math
import random
import numpy as np

random.seed(5)


def roots(deg):
    d = max(deg)
    c = np.zeros(d + 1)
    for k in deg:
        c[d - k] += 1
    return np.roots(c) if d > 0 else []


ok = True
worst, worst_real, cnt = 0.0, 0.0, 0
for n in range(2, 8):
    pairs = list(itertools.combinations(range(n), 2))
    E = len(pairs)
    masks = range(1 << E) if n <= 6 else (random.getrandbits(E) for _ in range(20000))
    for mask in masks:
        deg = [0] * n
        for i, (u, v) in enumerate(pairs):
            if mask >> i & 1:
                deg[u] += 1
                deg[v] += 1
        for z in roots(deg):
            if abs(z) < 1e-12:
                continue
            cnt += 1
            worst = max(worst, abs(z) / (n - 1) ** (abs(np.angle(z)) / math.pi))
            if abs(z.imag) < 1e-9:
                worst_real = max(worst_real, abs(z) / (n - 1))
print(f"{cnt} nonzero roots (all graphs on 2..6 vertices, 20000 random graphs on 7): max |z|/(n-1)^(|arg z|/pi) = {worst:.12f}; "
      f"real roots max |z|/(n-1) = {worst_real:.12f}")
ok &= worst <= 1 + 1e-9 and worst_real <= 1 + 1e-9
for n in range(3, 12):
    z = (n - 1) ** (1 / (n - 2)) * np.exp(1j * math.pi / (n - 2))
    res = abs(z ** (n - 1) + (n - 1) * z)
    ratio = abs(z) / (n - 1) ** (abs(np.angle(z)) / math.pi)
    print(f"star K_1,{n - 1}: |D(z)| = {res:.2e} at z = (n-1)^(1/(n-2)) e^(i pi/(n-2)); ratio = {ratio:.12f}")
    ok &= res < 1e-9 and abs(ratio - 1) < 1e-9
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
