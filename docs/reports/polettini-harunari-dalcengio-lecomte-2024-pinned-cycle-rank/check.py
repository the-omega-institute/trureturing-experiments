#!/usr/bin/env python3
"""Polettini-Harunari-Dal Cengio-Lecomte (arXiv:2407.16093v1; Lett. Math. Phys. 2025), conjecture (Introduction and
Conclusions): the 3^n spanning-tree vectors tau_x (one per root x) of the deletion-constriction through n pinned
reversible edges span an (n+1)-dimensional space.

tau_x(s), s in {0,+,-}^n: sum over spanning trees rooted at x (every edge directed toward x) of the product of rates,
restricted to trees containing neither direction of pinned edge i (s_i = 0), its forward direction only (+), or its
backward direction only (-).

(1) K4, all 12 rates 1, pins 1->2, 2->3, 3->1 (deleting them leaves the bidirected star at 4): exact integer rows
    as in the issue, R1 + R2 + R3 = 3 R4, rank 3 != 4.
(2) Random rational rates on K4, K5, K6 with triangle pins: rank 3 (exact rank over Q).
(3) The characterisation dim = 1 + rank(incidence of the pins): random graphs and random pin sets (forests and
    non-forests), exact ranks. Forest pins give n + 1.
Exit 0 and RESULT PASS.
"""
import itertools, random, sys
from fractions import Fraction as F
random.seed(20261019)

def trees(N, root, rate):
    """All spanning trees rooted at root: parent maps on the other vertices along positive-rate arcs, acyclic."""
    others = [v for v in range(N) if v != root]
    for par in itertools.product(range(N), repeat=len(others)):
        p = dict(zip(others, par))
        if any(p[v] == v or rate[v][p[v]] == 0 for v in others): continue
        ok = True
        for v in others:
            seen = set(); u = v
            while u != root:
                if u in seen: ok = False; break
                seen.add(u); u = p[u]
            if not ok: break
        if ok: yield p

def tau(N, rate, pins):
    rows = []
    for x in range(N):
        vec = {s: F(0) for s in itertools.product('0+-', repeat=len(pins))}
        for p in trees(N, x, rate):
            s = []
            for (u, v) in pins:
                fwd = p.get(u) == v; bwd = p.get(v) == u
                s.append('+' if fwd else '-' if bwd else '0')
            w = F(1)
            for a, b in p.items(): w *= rate[a][b]
            vec[tuple(s)] += w
        rows.append([vec[s] for s in itertools.product('0+-', repeat=len(pins))])
    return rows

def rank(M):
    M = [r[:] for r in M]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def inc_rank(N, pins):
    B = [[(1 if j == v else -1 if j == u else 0) for j in range(N)] for (u, v) in pins]
    return rank([[F(x) for x in row] for row in B])

ok = True
# (1)
N = 4; rate = [[F(0) if i == j else F(1) for j in range(N)] for i in range(N)]
pins = [(0, 1), (1, 2), (2, 0)]
R = tau(N, rate, pins)
expect = {
 0: [1,2,0,1,3,0,1,0,0,0,0,0,0,0,0,0,0,0,2,3,0,0,0,0,3,0,0],
 1: [1,1,1,0,0,0,2,0,3,2,3,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0],
 2: [1,0,2,2,0,3,0,0,0,1,0,0,3,0,0,0,0,0,1,0,3,0,0,0,0,0,0],
 3: [1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,0,0]}
ok &= all([int(v) for v in R[x]] == expect[x] for x in range(4))
ok &= all(R[0][k] + R[1][k] + R[2][k] == 3 * R[3][k] for k in range(27))
ok &= rank(R) == 3
print(f"K4 triangle pins: rank {rank(R)} (n + 1 = 4); R1+R2+R3 = 3 R4: {all(R[0][k]+R[1][k]+R[2][k]==3*R[3][k] for k in range(27))}")
# (2)
for N in [4, 5, 6]:
    for trial in range(2):
        rate = [[F(0) if i == j else F(random.randint(1, 9), random.randint(1, 9)) for j in range(N)] for i in range(N)]
        rk = rank(tau(N, rate, [(0, 1), (1, 2), (2, 0)]))
        ok &= rk == 3; print(f"K{N} random rates, triangle pins: rank {rk}")
# (3)
for trial in range(12):
    N = random.choice([4, 5, 6]); n = random.choice([2, 3])
    rate = [[F(0) if i == j else F(random.randint(1, 9), random.randint(1, 9)) for j in range(N)] for i in range(N)]
    pairs = random.sample([(i, j) for i in range(N) for j in range(i + 1, N)], n)
    rk = rank(tau(N, rate, pairs)); pred = 1 + inc_rank(N, pairs)
    ok &= rk == pred
    print(f"K{N} pins {pairs}: rank {rk}, 1 + rank(B) = {pred}, n + 1 = {n + 1}")
# (4) a pinned 4-cycle on K5 (n = 4): predicted 1 + rank(B) = 4, not n + 1 = 5
rate = [[F(0) if i == j else F(random.randint(1, 9), random.randint(1, 9)) for j in range(5)] for i in range(5)]
cyc = [(0, 1), (1, 2), (2, 3), (3, 0)]
rk = rank(tau(5, rate, cyc)); ok &= rk == 4 == 1 + inc_rank(5, cyc)
print(f"K5 random rates, pinned 4-cycle: rank {rk}, 1 + rank(B) = {1 + inc_rank(5, cyc)}, n + 1 = 5")
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
