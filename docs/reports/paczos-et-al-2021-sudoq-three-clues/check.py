#!/usr/bin/env python3
"""Checks for Paczos-Wierzbinski-Rajchel-Mieldzioc-Burchardt-Zyczkowski, arXiv:2106.02967v2,
Conjecture (Section III.A): no 4x4 SudoQ grid with 3 clues is uniquely solvable.
Usage: python3 check.py   (numpy).  Exits 0 and prints `bad= 0` when every check passes."""
import itertools, sys
import numpy as np

bad = 0
def ok(cond, msg):
    global bad
    if not cond:
        bad += 1
        print("FAIL", msg)

rng = np.random.default_rng(20261011)
CELLS = [(i, j) for i in range(4) for j in range(4)]
def block(i, j):
    return (i // 2, j // 2)
UNITS = ([[(i, j) for j in range(4)] for i in range(4)] +
         [[(i, j) for i in range(4)] for j in range(4)] +
         [[(2 * bi + a, 2 * bj + b) for a in range(2) for b in range(2)] for bi in range(2) for bj in range(2)])
ok(len(UNITS) == 12, "twelve units")

# 1. Sudoku symmetries act transitively on transversals (one cell per row, column and block).
trans = []
for perm in itertools.permutations(range(4)):
    cells = [(i, perm[i]) for i in range(4)]
    if len({block(*c) for c in cells}) == 4:
        trans.append(frozenset(cells))
gens = []
def rowperm(p):
    return lambda c: (p[c[0]], c[1])
def colperm(p):
    return lambda c: (c[0], p[c[1]])
for p in ([1, 0, 2, 3], [0, 1, 3, 2], [2, 3, 0, 1]):
    gens.append(rowperm(p)); gens.append(colperm(p))
gens.append(lambda c: (c[1], c[0]))
canon = frozenset([(0, 0), (1, 2), (2, 1), (3, 3)])
ok(canon in trans, "canonical pattern is a transversal")
orbit = {canon}; frontier = [canon]
while frontier:
    s = frontier.pop()
    for g in gens:
        t = frozenset(g(c) for c in s)
        if t not in orbit:
            orbit.add(t); frontier.append(t)
ok(orbit == set(trans), f"transitive on transversals: orbit {len(orbit)} of {len(trans)}")
# the generators preserve the unit structure
unitsets = {frozenset(u) for u in UNITS}
for g in gens:
    ok({frozenset(g(c) for c in u) for u in unitsets} == unitsets, "generator preserves units")

# 2. Forced equalities: two units sharing three rays force equality of the fourth.
sym = {c: f"x{c[0]}{c[1]}" for c in CELLS}
for c in canon:
    sym[c] = "w"
parent = {}
def find(x):
    while parent.get(x, x) != x:
        x = parent[x]
    return x
def union(x, y):
    parent[find(x)] = find(y)
changed = True
while changed:
    changed = False
    for u, v in itertools.combinations(UNITS, 2):
        su = {find(sym[c]) for c in u}; sv = {find(sym[c]) for c in v}
        if len(su) == 4 and len(sv) == 4 and len(su & sv) == 3:
            (x,) = su - sv; (y,) = sv - su
            if find(x) != find(y):
                union(x, y); changed = True
form = [[find(sym[(i, j)]) for j in range(4)] for i in range(4)]
names = {}
for row in form:
    for s in row:
        names.setdefault(s, len(names))
pattern = [[names[s] for s in row] for row in form]
# expected: [[w,a,b,c],[d,e,w,a],[g,w,h,b],[j,d,g,w]]
exp = [["w", "a", "b", "c"], ["d", "e", "w", "a"], ["g", "w", "h", "b"], ["j", "d", "g", "w"]]
en = {}
for row in exp:
    for s in row:
        en.setdefault(s, len(en))
ok(pattern == [[en[s] for s in row] for row in exp], f"forced form {pattern}")

# 3. The Q0 family and its four trades.
def haar(n):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / abs(np.diag(r)))
def is_onb(vs):
    M = np.array(vs).T
    return np.allclose(M.conj().T @ M, np.eye(4), atol=1e-9)
def valid(G):
    return all(is_onb([G[c] for c in u]) for u in UNITS)
def same_ray(x, y):
    return abs(abs(np.vdot(x, y)) - 1) < 1e-9
def q0(transpose=False):
    U = haar(4); w, a, b, c = U[:, 0], U[:, 1], U[:, 2], U[:, 3]
    V = haar(2); g = V[0, 0] * a + V[1, 0] * c; h = V[0, 1] * a + V[1, 1] * c
    rows = [[w, a, b, c], [b, c, w, a], [g, w, h, b], [h, b, g, w]]
    G = {(i, j): rows[i][j] for i in range(4) for j in range(4)}
    if transpose:
        G = {(i, j): G[(j, i)] for i in range(4) for j in range(4)}
    return G
RECTS = [[(r, c) for r in R for c in C] for R in ((0, 1), (2, 3)) for C in ((0, 2), (1, 3))]
for trial in range(300):
    tr = trial % 2 == 1
    G = q0(tr)
    ok(valid(G), "Q0 valid")
    rects = [[(c[1], c[0]) for c in r] for r in RECTS] if tr else RECTS
    cover = sorted(c for r in rects for c in r)
    ok(cover == sorted(CELLS), "rectangles partition")
    for r in rects:
        rays = []
        for c in r:
            if not any(same_ray(G[c], x) for x in rays):
                rays.append(G[c])
        ok(len(rays) == 2, "two rays per rectangle")
        x, y = rays
        H = dict(G)
        for c in r:
            H[c] = y if same_ray(G[c], x) else x
        ok(valid(H), "trade preserves constraints")
        ok(all(not same_ray(H[c], G[c]) for c in r), "trade changes every cell of the rectangle")

# 4. The constructive second completion on random grids with 3 random clues.
def sudoku_grids():
    base = [[0, 1, 2, 3], [2, 3, 0, 1], [1, 0, 3, 2], [3, 2, 1, 0]]
    while True:
        U = haar(4)
        p = rng.permutation(4)
        yield {(i, j): U[:, p[base[i][j]]] for i in range(4) for j in range(4)}
        yield q0(bool(rng.integers(2)))
def second_completion(G, clues):
    Svec = np.array([G[c] for c in clues]).T
    u, s, vh = np.linalg.svd(Svec)
    rank = int((s > 1e-9).sum())
    PS = u[:, :rank] @ u[:, :rank].conj().T
    T = 2 * PS - np.eye(4)
    H = {c: T @ G[c] for c in CELLS}
    if any(not same_ray(H[c], G[c]) for c in CELLS):
        return H
    K = [c for c in CELLS if np.linalg.norm(PS @ G[c]) < 1e-9]
    if 4 - rank >= 2:
        row0 = [G[(0, j)] for j in range(4) if any((0, j) == k for k in K)]
        x, y = row0[0], row0[1]
        Uk = np.eye(4) - np.outer(x, x.conj()) - np.outer(y, y.conj()) + np.outer(x, y.conj()) + np.outer(y, x.conj())
        return {c: Uk @ G[c] for c in CELLS}
    for r in RECTS + [[(c[1], c[0]) for c in r] for r in RECTS]:
        if any(c in clues for c in r):
            continue
        rays = []
        for c in r:
            if not any(same_ray(G[c], x) for x in rays):
                rays.append(G[c])
        if len(rays) != 2:
            continue
        x, y = rays
        H = dict(G)
        for c in r:
            H[c] = y if same_ray(G[c], x) else x
        if valid(H):
            return H
    return None
gen = sudoku_grids()
for trial in range(2000):
    G = next(gen)
    ok(valid(G), "random grid valid")
    clues = [CELLS[k] for k in rng.choice(16, size=3, replace=False)]
    H = second_completion(G, clues)
    ok(H is not None, "second completion found")
    if H is None:
        continue
    ok(valid(H), "second completion valid")
    ok(all(same_ray(H[c], G[c]) for c in clues), "clues kept")
    ok(any(not same_ray(H[c], G[c]) for c in CELLS), "completion differs")

print("transversals", len(trans))
print("bad=", bad)
sys.exit(1 if bad else 0)
