"""Juhnke–Poullot (arXiv:2504.20739v3) Problem 3.12 checks under the binding conventions of the preregistration.

Edges: [u, v] is an edge iff some c makes u, v the only maximizers of <x, c> over the vertex set
(exact rational LP via sympy). Monotone paths: directed paths in G_{P,c} from v_min to v_max.
Coherent paths: upper paths of the projections (<x,c>, <x,w>) for random w (sampled; lower bound on the set).
Pyr(P) = conv({0} u P x {1}), c_{k+1} = (c_k, 1 - min <x, c_k>).
"""
import itertools, random, sys
from fractions import Fraction as Fr
from collections import Counter
from sympy import symbols, Rational
from sympy.solvers.simplex import lpmax
import numpy as np

def dot(x, c): return sum(a * b for a, b in zip(x, c))

def edges(V):
    n, d = len(V), len(V[0])
    cs = symbols(f'c0:{d}')
    E = set()
    for i, j in itertools.combinations(range(n), 2):
        u, v = V[i], V[j]
        cons = [sum(cs[t] * (Rational(u[t]) - Rational(v[t])) for t in range(d)) <= 0,
                sum(cs[t] * (Rational(v[t]) - Rational(u[t])) for t in range(d)) <= 0]
        for k in range(n):
            if k in (i, j): continue
            cons.append(sum(cs[t] * (Rational(u[t]) - Rational(V[k][t])) for t in range(d)) >= 1)
        cons += [cs[t] <= 1000 for t in range(d)] + [cs[t] >= -1000 for t in range(d)]
        try:
            lpmax(Rational(0), cons); E.add((i, j))
        except Exception:
            pass
    return E

def monotone_counts(V, c, E):
    h = [dot(v, c) for v in V]
    for i, j in E: assert h[i] != h[j], 'c not generic'
    adj = {i: [] for i in range(len(V))}
    for i, j in E:
        (a, b) = (i, j) if h[i] < h[j] else (j, i)
        adj[a].append(b)
    src = min(range(len(V)), key=lambda i: h[i]); snk = max(range(len(V)), key=lambda i: h[i])
    cnt = Counter()
    def go(x, l):
        if x == snk: cnt[l] += 1; return
        for y in adj[x]: go(y, l + 1)
    go(src, 0)
    return [cnt.get(l, 0) for l in range(max(cnt) + 1)]

rng = np.random.default_rng(7)
def coherent_counts(V, c, samples):
    A = np.array([[float(t) for t in v] for v in V]); h = A @ np.array([float(t) for t in c])
    order = np.argsort(h); paths = set()
    for _ in range(samples // 2000):
        Y = A @ rng.standard_normal((A.shape[1], 2000))
        for s in range(Y.shape[1]):
            hull = []
            for i in order:
                p = (h[i], Y[i, s], i)
                while len(hull) >= 2:
                    (x1, y1, _), (x2, y2, _) = hull[-2], hull[-1]
                    if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) >= 0: hull.pop()
                    else: break
                hull.append(p)
            paths.add(tuple(q[2] for q in hull))
    cnt = Counter(len(p) - 1 for p in paths)
    return [cnt.get(l, 0) for l in range(max(cnt) + 1)]

def pyr(V, c):
    m = min(dot(v, c) for v in V)
    return [tuple([0] * (len(c) + 1))] + [tuple(v) + (1,) for v in V], tuple(c) + (1 - m,)

def unimodal(a):
    a = list(a); i = 0
    while i + 1 < len(a) and a[i + 1] >= a[i]: i += 1
    while i + 1 < len(a) and a[i + 1] <= a[i]: i += 1
    return i == len(a) - 1

def mul_one_plus_x(a, m):
    out = list(a)
    for _ in range(m): out = [x + y for x, y in zip(out + [0], [0] + out)]
    return out

cube = [tuple(Fr(t) for t in v) for v in itertools.product([0, 1], repeat=3)]
lop3 = [tuple(map(Fr, v)) for v in [(0,0,0),(0,1,0),(0,0,1),(0,Fr(1,3),1),(1,0,0),(4,1,0),(2,0,1),(3,Fr(1,3),1)]]
ok = True
for name, V0, c0, expect in [('cube3', cube, (Fr(1), Fr(2), Fr(4)), [0, 0, 0, 6]), ('Lop3', lop3, (Fr(1),) * 3, [0, 0, 2, 0, 4])]:
    V, c = V0, c0; base = None
    for k in range(0, 4):
        if k: V, c = pyr(V, c)
        E = edges(V)
        mono = monotone_counts(V, c, E); coh = coherent_counts(V, c, 200000)
        print(f'{name} k={k}: |V|={len(V)} |E|={len(E)} mono={mono} coh={coh} unimodal=({unimodal(mono)},{unimodal(coh)})')
        if k == 0:
            print(f'  fidelity (Example 3.5 / Example 4.1): mono == coh == {expect}: {mono == expect and coh == expect}')
            ok &= mono == expect and coh == expect
        if k == 1: base = (mono, coh)
        if k >= 2:
            pm = mul_one_plus_x(base[0], k - 1); pc = mul_one_plus_x(base[1], k - 1)
            print(f'  product structure X(1+X)^(k-1) Q: mono {mono == pm}, coh {coh == pc}')
            ok &= mono == pm and coh == pc
print('ALL_CHECKS_OK' if ok else 'CHECK_FAILED'); sys.exit(0 if ok else 1)
