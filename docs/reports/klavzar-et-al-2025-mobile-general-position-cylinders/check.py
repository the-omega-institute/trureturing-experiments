"""Checks for Klavzar, Krishnakumar, Kuziak, Shallcross, Tuite and Yero, arXiv:2505.00535v2 (DAM 379 (2026) 768-780),
the conjecture after Proposition 3.7: Mob_gp(C_7 x P_s) = Mob_gp(C_10 x P_s) = 4 for s >= 5.
Distances come from breadth-first search on the Cartesian product graph (no formula assumed).
(1) Five-point angular rigidity: the words (x_0..x_4) in Z_r with every ordered triple strictly off a shortest cycle path
    are exactly b + e*a*i (r = 3a + 1): 14 words for r = 7, 20 for r = 10.
(2) For s = 2..6 no six vertices of C_r x P_s are in general position (exhaustive, with pruning).
(3) For s = 2..6 every five-vertex general-position configuration has no legal move that changes a cycle coordinate,
    and every legal move keeps the set of occupied cycle coordinates (so some cycle column is never visited).
(4) The explicit four-robot traversal (two-row cycle sweep plus row shifts) is legal at every step and visits every
    vertex for s = 2..8.  Prints ALL_OK."""
import itertools
from collections import deque
def graph(r, s):
    V = [(x, y) for x in range(r) for y in range(s)]
    adj = {v: [] for v in V}
    for (x, y) in V:
        adj[(x, y)] += [((x + 1) % r, y), ((x - 1) % r, y)] + ([(x, y + 1)] if y + 1 < s else []) + ([(x, y - 1)] if y > 0 else [])
    dist = {}
    for v in V:
        d = {v: 0}; q = deque([v])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in d: d[w] = d[u] + 1; q.append(w)
        dist[v] = d
    return V, adj, dist
def between(dist, u, v, w):  # v on a shortest u-w path
    return dist[u][v] + dist[v][w] == dist[u][w]
def gp(dist, S):
    S = list(S)
    for u, v, w in itertools.permutations(S, 3):
        if u < w and between(dist, u, v, w): return False
    return True
ok = True
# (1)
for a in (2, 3):
    r = 3 * a + 1
    dc = lambda x, z: min((x - z) % r, (z - x) % r)
    words = [w for w in itertools.product(range(r), repeat=5)
             if all(dc(w[i], w[j]) + dc(w[j], w[k]) > dc(w[i], w[k]) for i, j, k in itertools.combinations(range(5), 3))]
    ap = {tuple((b + e * a * i) % r for i in range(5)) for b in range(r) for e in (1, -1)}
    ok &= set(words) == ap and len(words) == {7: 14, 10: 20}[r]
print("(1) five-point angular rigidity:", ok)
# (2), (3)
ok2 = True; counts = {}
for r in (7, 10):
    for s in range(2, 7):
        V, adj, dist = graph(r, s)
        idx = {v: i for i, v in enumerate(V)}
        # extend general-position sets incrementally (sets of size k, sorted)
        level = [(v,) for v in V]
        five = []
        for k in range(2, 7):
            nxt = []
            for S in level:
                for v in V:
                    if idx[v] <= idx[S[-1]]: continue
                    T = S + (v,)
                    if all(not between(dist, p, q, t) for p, q, t in itertools.permutations(T, 3) if v in (p, q, t) and p < t):
                        nxt.append(T)
            level = nxt
            if k == 5: five = list(level)
            if k == 6: ok2 &= len(level) == 0
        counts[(r, s)] = len(five)
        for S in five:
            Sset = set(S)
            for u in S:
                for v in adj[u]:
                    if v in Sset: continue
                    T = (Sset - {u}) | {v}
                    if gp(dist, T):
                        ok2 &= v[0] == u[0]      # only vertical legal moves
print("(2) no six vertices in general position, (3) five robots move only vertically, s = 2..6:", ok2, "five-sets:", counts)
# (4) explicit traversal
def traversal(r, s):
    a = (r - 1) // 3
    L0, H0 = ({0, 4}, {1, 5}) if r == 7 else ({0, 5}, {1, 6})
    seq7 = [("H", 1, 2), ("L", 0, 1), ("H", 5, 6), ("L", 4, 5)]
    seq10 = [("H", 1, 2), ("H", 6, 7), ("L", 0, 1), ("L", 5, 6)]
    seq = seq7 if r == 7 else seq10
    conf = {(x, 0) for x in L0} | {(x, 1) for x in H0}
    states = [frozenset(conf)]; visited = set(conf)
    def move(u, v):
        nonlocal conf
        assert u in conf and v not in conf
        conf = (conf - {u}) | {v}; states.append(frozenset(conf)); visited.add(v)
    for j in range(s - 1):
        if j > 0:
            ups = sorted([p for p in conf if p[1] == j], key=lambda p: p[0])
            lows = sorted([p for p in conf if p[1] == j - 1], key=lambda p: p[0])
            for p in ups: move(p, (p[0], j + 1))
            for p in lows: move(p, (p[0], j))
        for t in range(r):
            for row, x0, x1 in seq:
                y = j + 1 if row == "H" else j
                move(((x0 + t) % r, y), ((x1 + t) % r, y))
    return states, visited
ok3 = True
for r in (7, 10):
    for s in range(2, 9):
        V, adj, dist = graph(r, s)
        states, visited = traversal(r, s)
        ok3 &= all(gp(dist, S) and len(S) == 4 for S in states) and visited == set(V)
        ok3 &= all(len(A ^ B) == 2 and any(b in adj[a] for a in A - B for b in B - A) for A, B in zip(states, states[1:]))
print("(4) four-robot traversal legal and covering for s = 2..8:", ok3)
print("ALL_OK" if ok and ok2 and ok3 else "FAIL")
