#!/usr/bin/env python3
"""Exact checks for the subcubic partial-Petrial counterexample.

All polynomial coefficients and GF(2) ranks are computed with integers.  The
47-vertex computation uses the rooted two-state recurrence derived from the
matrix-rank definition; small trees are independently enumerated by binary
Gaussian elimination.
"""
from functools import lru_cache
import random


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return trim(c)


def psub(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] -= x
    return trim(c)


def pscale(a, k):
    return trim([k * x for x in a])


def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def pshift(a, k):
    return [0] * k + list(a)


def rank_gf2(rows, n):
    """Rank of n x n binary rows represented by integer bit masks."""
    rows = list(rows)
    rank = 0
    for col in range(n):
        pivot = next((i for i in range(rank, n) if (rows[i] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(n):
            if i != rank and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def direct_polynomial(n, edges):
    """P_G^times from the literal diagonal-perturbation rank sum."""
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    counts = [0] * (n + 1)
    for diagonal in range(1 << n):
        rows = [adj[i] ^ (((diagonal >> i) & 1) << i) for i in range(n)]
        counts[rank_gf2(rows, n)] += 1
    return trim(counts)


def rooted_polynomial(n, edges, root=0):
    """Two-state rooted recurrence for a tree, with the matrix definition below it."""
    children = [[] for _ in range(n)]
    seen = {root}

    def orient(v):
        for a, b in edges:
            if a == v and b not in seen:
                seen.add(b)
                children[v].append(b)
                orient(b)
            elif b == v and a not in seen:
                seen.add(a)
                children[v].append(a)
                orient(a)

    orient(root)

    @lru_cache(None)
    def states(v):
        x, y = [1], [1]
        for child in children[v]:
            a, b = states(child)
            x = pmul(x, padd(a, b))
            y = pmul(y, b)
        # A: some kernel vector is nonzero at the root; B: every kernel
        # vector vanishes there.  The parent matrix has diagonal 0 or 1.
        return y, padd(pshift(y, 1), pshift(pscale(psub(x, y), 2), 2))

    a, b = states(root)
    return padd(a, b)


BASE_EDGES = [
    (0,1),(1,2),(0,3),(2,4),(3,5),(3,6),(0,7),(7,8),(2,9),(4,10),
    (5,11),(10,12),(6,13),(8,14),(10,15),(11,16),(6,17),(5,18),(13,19),
    (4,20),(13,21),(11,22),(7,23),(8,24),(21,25),(17,26),(1,27),(26,28),
    (16,29),(19,30),(12,31),(20,32),(23,33),(32,34),(24,35),(24,36),(12,37),
    (27,38),(20,39),(33,40),(27,41),(41,42),(41,43),(23,44),(32,45),(33,46),
]
EXPECTED_P = [
    10240,1627136,260295680,6854209024,86056429568,644755355008,
    3154826928992,10490698794208,23971970245256,37152671631956,
    37379380458544,22053301655367,5796710714349,
]
EXPECTED_Q = [
    173056,102174208,2890419712,37693663872,291483305024,1467484077824,
    5006704361856,11705320102888,18511478712388,18957977106970,
    11360447947371,3027162132495,
]


def make_tail(r):
    edges = list(BASE_EDGES)
    if r:
        prev = 22
        for v in range(47, 47 + r):
            edges.append((prev, v))
            prev = v
    return 47 + r, edges


def prufer_tree(n, rng):
    if n == 1:
        return []
    seq = [rng.randrange(n) for _ in range(n - 2)]
    degree = [1] * n
    for x in seq:
        degree[x] += 1
    edges = []
    for x in seq:
        leaf = next(i for i, d in enumerate(degree) if d == 1)
        edges.append((leaf, x))
        degree[leaf] -= 1
        degree[x] -= 1
    left = [i for i, d in enumerate(degree) if d == 1]
    edges.append((left[0], left[1]))
    return edges


def jacobsthal(k):
    vals = [0, 1]
    for _ in range(k - 1):
        vals.append(vals[-1] + 2 * vals[-2])
    return vals[k]


def main():
    assert len(BASE_EDGES) == 46
    degree = [0] * 47
    for u, v in BASE_EDGES:
        degree[u] += 1
        degree[v] += 1
    assert sorted({*degree}) == [1, 2, 3]
    assert degree.count(1) == 22 and degree.count(2) == 5 and degree.count(3) == 20
    print("base tree: 47 vertices, 46 edges, degree counts (1,2,3)=(22,5,20)")

    # Independent literal rank enumeration on 100 small trees.
    rng = random.Random(2611)
    checked = 0
    for n in range(1, 11):
        for _ in range(10):
            edges = prufer_tree(n, rng)
            assert rooted_polynomial(n, edges) == direct_polynomial(n, edges)
            checked += 1
    print(f"small-tree rank-definition cross-checks: {checked}")

    # Rooted states evaluate the two 47/46 vertex polynomials.
    p = rooted_polynomial(47, BASE_EDGES)
    q_edges = [e for e in BASE_EDGES if 22 not in e]
    # Removing leaf 22 leaves 46 vertices; relabel the remaining graph.
    rem = [v for v in range(47) if v != 22]
    ren = {v: i for i, v in enumerate(rem)}
    q = rooted_polynomial(46, [(ren[u], ren[v]) for u, v in q_edges])
    assert p[35:48] == EXPECTED_P and all(x == 0 for x in p[:35])
    assert q[35:47] == EXPECTED_Q and all(x == 0 for x in q[:35])
    assert sum(p) == 2 ** 47 and sum(q) == 2 ** 46
    deficit = p[36] * p[36] - p[35] * p[37]
    assert deficit == -17856200704
    print("rooted-state P coefficients 35..47:", p[35:48])
    print("rooted-state Q coefficients 35..46:", q[35:47])
    print("sum(P), sum(Q), c36^2-c35*c37:", sum(p), sum(q), deficit)

    # Verify the tail identity on actual rooted trees for initial r, then
    # check its exact Jacobsthal inequality for 0 <= r <= 100.
    for r in range(7):
        n, edges = make_tail(r)
        actual = rooted_polynomial(n, edges)
        a, b = jacobsthal(r + 1), jacobsthal(r)
        formula = pshift(padd(pscale(p, a), pshift(pscale(q, 2 * b), 1)), r)
        assert actual == formula
    for r in range(101):
        a, b = jacobsthal(r + 1), jacobsthal(r)
        formula = pshift(padd(pscale(p, a), pshift(pscale(q, 2 * b), 1)), r)
        assert sum(formula) == 2 ** (47 + r)
        lhs = (1024 * (1589 * a + 338 * b)) ** 2
        rhs = (1024 * (10 * a)) * (1024 * (254195 * a + 199559 * b))
        assert lhs < rhs
        assert 0 <= b <= a
    print("tail identity: rooted-tree checks r=0..6; exact family inequalities r=0..100")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
