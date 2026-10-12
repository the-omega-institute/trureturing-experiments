#!/usr/bin/env python3
"""Exact check of the kappa^2 entropy-ordering question of Atashpendar-Mestel-Roscoe-Ryan (arXiv:1807.11609v1, Sec. V).
For binary x of length m and n > m: omega_x(y) = number of index sets pi of size m with y_pi = x;
H_n(x) = -sum_y p log2 p with p = omega_x(y) / (C(n,m) 2^(n-m)); kappa2(x) = sum_{r,s} C(r+s-2,r-1) C(2m-r-s,m-r) [x_r = x_s].
Reports every strict inversion kappa2(u) > kappa2(v) with H_n(u) > H_n(v) (entropies compared exactly through
sum_y w log w, since H = log2 mu - (1/mu) sum w log2 w with mu common to all x of the same (n, m))."""
import itertools, math, sys
from fractions import Fraction

def emb(x, y):
    m = len(x); D = [1] + [0] * m
    for c in y:
        for i in range(m, 0, -1):
            if x[i - 1] == c:
                D[i] += D[i - 1]
    return D[m]

def kappa2(x):
    m = len(x)
    return sum(math.comb(r + s - 2, r - 1) * math.comb(2 * m - r - s, m - r)
               for r in range(1, m + 1) for s in range(1, m + 1) if x[r - 1] == x[s - 1])

def weights(x, n):
    return [emb(x, y) for y in itertools.product((0, 1), repeat=n)]

def S(ws):
    # sum w log w as an exact multiplicative object: prod w^w (compare H via these, larger product = smaller H)
    p = 1
    for w in ws:
        if w > 1:
            p *= w ** w
    return p

def run(n, m):
    mu = math.comb(n, m) * 2 ** (n - m)
    rows = []
    for x in itertools.product((0, 1), repeat=m):
        ws = weights(x, n)
        assert sum(ws) == mu
        rows.append((''.join(map(str, x)), kappa2(x), S(ws), ws))
    inv = []
    for (a, ka, sa, _), (b, kb, sb, _) in itertools.permutations(rows, 2):
        # kappa(a) > kappa(b) predicts H(a) <= H(b), i.e. prod_a w^w >= prod_b w^w
        if ka > kb and sa < sb:
            inv.append((a, ka, b, kb))
    return rows, inv

if __name__ == "__main__":
    bad_total = 0
    for (n, m) in [(6, 5), (8, 5), (7, 5), (9, 5), (6, 4), (8, 6)]:
        rows, inv = run(n, m)
        H = {r[0]: math.log2(math.comb(n, m) * 2 ** (n - m)) - sum(w * math.log2(w) for w in r[3] if w > 0) / (math.comb(n, m) * 2 ** (n - m)) for r in rows}
        k = {r[0]: r[1] for r in rows}
        print(f"n={n} m={m} strict_inversions={len(inv)}")
        for a, ka, b, kb in inv[:4]:
            print(f"   kappa2({a})={ka} > kappa2({b})={kb} but H({a})={H[a]:.6f} > H({b})={H[b]:.6f}")
        if (n, m) in [(6, 5), (8, 5)]:
            print(f"   kappa2(01110)={k['01110']} kappa2(00100)={k['00100']} H(01110)={H['01110']:.6f} H(00100)={H['00100']:.6f} kappa2(11111)={k['11111']} H(11111)={H['11111']:.4f} H(01010)={H['01010']:.4f}")
    # family x_r = 0 1^(2r-1) 0, y_r = 0^r 1 0^r, m = 2r+1, n = 2r+2
    for r in range(2, 6):
        m, n = 2 * r + 1, 2 * r + 2
        x = (0,) + (1,) * (2 * r - 1) + (0,); y = (0,) * r + (1,) + (0,) * r
        kx, ky = kappa2(x), kappa2(y); sx, sy = S(weights(x, n)), S(weights(y, n))
        ok = ky > kx and sy < sx
        print(f"family r={r} m={m} n={n} kappa2(x)={kx} kappa2(y)={ky} H(y)>H(x):{sy < sx} inversion:{ok}")
