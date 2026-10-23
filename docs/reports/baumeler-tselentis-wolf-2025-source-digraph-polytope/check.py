#!/usr/bin/env python3
"""Checks for A. Baumeler, E.-E. Tselentis and S. Wolf, "Graphical Tests of Causality", arXiv:2511.22552v1, Section 5:
the identity Q~_n^src = Q_n^src after Definition 13.

E = {(s, r) : s != r} on [n]; alpha(D) is the 0/1 adjacency vector; Q = conv{alpha(D) : D has a vertex of in-degree 0}
(Definition 8, Theorem 3); Q~ is cut out by x_sr >= 0, x_sr <= 1 (n >= 3) and sum_r x_{f(r) r} <= n - 1 for every
fixed-point-free f (Definitions 9 and 13).
Fidelity: Theorem 4 (every kefalopoda inequality is facet-defining) for n = 3, 4 by the rank of the tight source digraphs,
with a count of the f for which the proof's set D_{i,j} has repeated members; the two non-isomorphic kefalopoda digraphs for n = 3; Q~ = Q for n = 3 by vertex enumeration.
Proof checks (exact rationals): Q is inside Q~ for n <= 4 exhaustively; for random points of Q~ (n = 2..8, with
boundary cases), the column-maximum form of the kefalopoda inequalities (n <= 6, all f) and the explicit convex
decomposition into fixed-source boxes. Exits 1 on failure.
"""
import itertools, random, sys
from fractions import Fraction as Fr
import numpy as np

ok = True
def arcs(n): return [(s, r) for s in range(n) for r in range(n) if s != r]
def fpf(n): return [f for f in itertools.product(range(n), repeat=n) if all(f[i] != i for i in range(n))]
def has_source(n, A): return any(all((s, i) not in A for s in range(n) if s != i) for i in range(n))
# fidelity: Theorem 4 (each kefalopoda inequality is facet-defining) by direct rank computation, and the proof's set D_ij
for n in (3, 4):
    E = arcs(n); d = len(E)
    srcs = [A for k in range(d + 1) for A in map(set, itertools.combinations(E, k)) if has_source(n, A)]
    facet = True; dup = 0
    for f in fpf(n):
        kf = {(f[i], i) for i in range(n)}
        tight = [[1 if e in A else 0 for e in E] + [1] for A in srcs if sum(1 for e in kf if e in A) == n - 1]
        facet &= np.linalg.matrix_rank(np.array(tight, float)) == d
        Dij = {frozenset(({(i, j)} | kf) - {(f[i], i)}) for (i, j) in E}
        dup += len(Dij) < d
    print(f"fidelity Theorem 4, n = {n}: every kefalopoda inequality is facet-defining (tight source digraphs span a hyperplane) for all {len(fpf(n))} f: {facet}; f for which the proof's digraphs D_ij are not {d} distinct: {dup}")
    ok &= facet
def canon(n, f):
    return min(tuple(sorted((p.index(f[i]), p.index(i)) for i in range(n))) for p in itertools.permutations(range(n)))
iso3 = len({canon(3, f) for f in fpf(3)})
print(f"fidelity: non-isomorphic kefalopoda digraphs for n = 3: {iso3}")
ok &= iso3 == 2
# fidelity: vertex enumeration of Q~ for n = 3
n = 3; E = arcs(n); d = len(E)
rows = [([-1 if e == a else 0 for e in E], 0) for a in E] + [([1 if e == a else 0 for e in E], 1) for a in E]
rows += [([1 if e in {(f[i], i) for i in range(n)} else 0 for e in E], n - 1) for f in fpf(n)]
Aq = np.array([r[0] for r in rows], float); bq = np.array([r[1] for r in rows], float)
verts = set()
for T in itertools.combinations(range(len(rows)), d):
    M = Aq[list(T)]
    if abs(np.linalg.det(M)) < 1e-9: continue
    x = np.linalg.solve(M, bq[list(T)])
    if np.all(Aq @ x <= bq + 1e-9): verts.add(tuple(np.round(x, 9)))
src_vecs = {tuple(float(e in A) for e in E) for k in range(d + 1) for A in map(set, itertools.combinations(E, k)) if has_source(n, A)}
v_ok = all(all(abs(c - round(c)) < 1e-9 for c in v) for v in verts) and {tuple(float(round(c)) for c in v) for v in verts} <= src_vecs
print(f"fidelity n = 3: Q~ has {len(verts)} vertices, all 0/1 adjacency vectors of source digraphs (so Q~ = Q): {v_ok}")
ok &= v_ok
# Q inside Q~, n <= 4 exhaustive
for n in (2, 3, 4):
    E = arcs(n); F = fpf(n); bad = 0
    for k in range(len(E) + 1):
        for A in map(set, itertools.combinations(E, k)):
            if not has_source(n, A): continue
            if any(sum(1 for i in range(n) if (f[i], i) in A) > n - 1 for f in F): bad += 1
    print(f"n = {n}: every source digraph satisfies every kefalopoda inequality: {bad == 0}")
    ok &= bad == 0
# random points of Q~ and the decomposition
rng = random.Random(20261023)
def random_point(n):
    E = arcs(n)
    budget = Fr(n - 1)
    m = {}
    order = list(range(n)); rng.shuffle(order)
    for r in order:
        mode = rng.random()
        hi = min(Fr(1), budget)
        m[r] = Fr(0) if mode < 0.15 else (hi if mode < 0.35 else Fr(rng.randint(0, 60), 60) * hi)
        budget -= m[r]
    x = {}
    for r in range(n):
        ss = [s for s in range(n) if s != r]; top = rng.choice(ss)
        for s in ss: x[(s, r)] = m[r] if s == top else Fr(rng.randint(0, 60), 60) * m[r]
    return x
dec_ok = True; lem_ok = True; cnt = 0
for n in range(2, 9):
    E = arcs(n); F = fpf(n) if n <= 6 else None
    for _ in range(300 if n <= 6 else 100):
        x = random_point(n); cnt += 1
        m = {r: max(x[(s, r)] for s in range(n) if s != r) for r in range(n)}
        assert sum(m.values()) <= n - 1 and all(0 <= v <= 1 for v in x.values())
        if F is not None:
            lem_ok &= max(sum(x[(f[r], r)] for r in range(n)) for f in F) == sum(m.values())
        c = {r: 1 - m[r] for r in range(n)}; C = sum(c.values()); qv = {r: c[r] / C for r in range(n)}
        y = {}
        for i in range(n):
            for (s, r) in E:
                y[(i, s, r)] = Fr(0) if r == i else (x[(s, r)] / (1 - qv[r]) if qv[r] < 1 else Fr(0))
        dec_ok &= C >= 1 and all(0 <= v <= 1 for v in y.values()) and all(y[(i, s, i)] == 0 for i in range(n) for s in range(n) if s != i)
        dec_ok &= all(sum(qv[i] * y[(i, s, r)] for i in range(n)) == x[(s, r)] for (s, r) in E) and sum(qv.values()) == 1 and all(v >= 0 for v in qv.values())
print(f"{cnt} random points of Q~ (n = 2..8): column-maximum form of the kefalopoda inequalities (n <= 6, all f) {lem_ok}; convex decomposition x = sum_i q_i y^(i) with y^(i) in the box with source i {dec_ok}")
ok &= lem_ok and dec_ok
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
