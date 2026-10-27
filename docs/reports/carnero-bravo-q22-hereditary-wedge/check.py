#!/usr/bin/env python3
"""Exact finite diagnostics for Carnero Bravo, arXiv:2505.06457v2, Q22.

Standard library only. Faces include the empty face (mask 0). Rational Betti
numbers are diagnostics, not a decision procedure for wedge homotopy type.
"""
from functools import lru_cache
from itertools import combinations
from math import gcd


def graph(n, edges):
    a = [0] * n
    for u, v in edges:
        assert u != v and 0 <= u < n and 0 <= v < n
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def all_graphs(n):
    pairs = list(combinations(range(n), 2))
    for e in range(1 << len(pairs)):
        yield graph(n, (p for i, p in enumerate(pairs) if e >> i & 1))


def bits(mask):
    return [i for i in range(mask.bit_length()) if mask >> i & 1]


@lru_cache(None)
def faces(a, allowed=None):
    if allowed is None:
        allowed = (1 << len(a)) - 1
    if not allowed:
        return frozenset({0})
    vbit = allowed & -allowed
    v = vbit.bit_length() - 1
    rest = allowed ^ vbit
    return faces(a, rest) | frozenset(vbit | f for f in faces(a, rest & ~a[v]))


def join(k, l):
    assert all(not (s & t) for s in k for t in l)
    return frozenset(s | t for s in k for t in l)


def rank_q(matrix):
    """Integer row elimination computing rank over Q; no float operations."""
    a = [row[:] for row in matrix]
    if not a:
        return 0
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        p = a[row][col]
        for i in range(row + 1, len(a)):
            q = a[i][col]
            if q:
                a[i] = [p * x - q * y for x, y in zip(a[i], a[row])]
                d = 0
                for z in a[i]:
                    d = gcd(d, z)
                if d:
                    a[i] = [z // d for z in a[i]]
        row += 1
        if row == len(a):
            break
    return row


def betti(k):
    """Ordinary rational Betti vector, omitting trailing zeroes."""
    levels = [[f for f in sorted(k) if f.bit_count() == q + 1]
              for q in range(max(f.bit_count() for f in k))]
    ranks = [0]
    for q in range(1, len(levels)):
        index = {s: i for i, s in enumerate(levels[q - 1])}
        mat = [[0] * len(levels[q]) for _ in levels[q - 1]]
        for j, s in enumerate(levels[q]):
            for i, v in enumerate(bits(s)):
                mat[index[s ^ (1 << v)]][j] = (-1) ** i
        ranks.append(rank_q(mat))
    ranks.append(0)
    b = [len(levels[q]) - ranks[q] - ranks[q + 1] for q in range(len(levels))]
    while b and b[-1] == 0:
        b.pop()
    return tuple(b)


def collapse(k, free, top):
    """Verify and perform an elementary simplicial collapse."""
    assert free and free in k and top in k
    assert free & top == free and top.bit_count() == free.bit_count() + 1
    assert {s for s in k if s & free == free} == {free, top}
    return k - {free, top}


def substitute(a, v, h):
    n, r = len(a), len(h)
    # Original v remains the first vertex of its replacement fiber.
    fiber = [v] + list(range(n, n + r - 1))
    old = [i for i in range(n) if i != v]
    edges = [(i, j) for i, j in combinations(old, 2) if a[i] >> j & 1]
    edges += [(fiber[i], fiber[j]) for i, j in combinations(range(r), 2) if h[i] >> j & 1]
    edges += [(i, w) for i in old for w in fiber if a[v] >> i & 1]
    return graph(n + r - 1, edges), fiber


def lex(a, h):
    n, r = len(a), len(h)
    return graph(n * r, ((i * r + x, j * r + y)
                         for i in range(n) for x in range(r)
                         for j in range(n) for y in range(r)
                         if i * r + x < j * r + y and
                         ((a[i] >> j & 1) or (i == j and (h[x] >> y & 1)))))


def source_fidelity():
    # Theorem 7, n=5=3*1+2: I(P5) is homotopy equivalent to S^1.
    p5 = graph(5, [(i, i + 1) for i in range(4)])
    k = faces(p5)
    k = collapse(k, (1 << 0) | (1 << 2), (1 << 0) | (1 << 2) | (1 << 4))
    k = collapse(k, 1 << 2, (1 << 2) | (1 << 4))
    square_edges = {frozenset({0, 3}), frozenset({3, 1}),
                    frozenset({1, 4}), frozenset({4, 0})}
    assert {frozenset(bits(s)) for s in k if s.bit_count() == 2} == square_edges
    assert max(s.bit_count() for s in k) == 2 and betti(k) == (1, 1)
    print('FIDELITY: Theorem 7, I(P5) collapses in 2 verified elementary steps to a 4-cycle; Betti=(1,1)')


def small_substitution_checks():
    cases = restrictions = link_empty = 0
    for n in range(1, 5):
        for a in all_graphs(n):
            for r in range(1, 4):
                for h in all_graphs(r):
                    for v in range(n):
                        s, fiber = substitute(a, v, h)
                        fiber_mask = sum(1 << w for w in fiber)
                        old_mask = ((1 << n) - 1) ^ (1 << v)
                        k = faces(a)
                        deletion = faces(a, old_mask)
                        link = faces(a, old_mask & ~a[v])
                        # Graph-link identity, checked without assuming homotopy type.
                        assert link == frozenset(f ^ (1 << v) for f in k if f >> v & 1)
                        x = frozenset(sum(1 << fiber[j] for j in bits(f)) for f in faces(h))
                        lx = join(link, x)
                        cone = join(link, frozenset({0, 1 << v}))
                        assert deletion | cone == k
                        assert k & lx == cone
                        assert faces(s) == deletion | lx
                        if link == frozenset({0}):
                            link_empty += 1
                        else:
                            # Ordinary reduced-homology consequence of the connected-component splitting.
                            b, c, d = betti(faces(s)), betti(k), betti(lx)
                            assert b[0] == c[0]
                            for q in range(1, max(len(b), len(c), len(d))):
                                at = lambda z: z[q] if q < len(z) else 0
                                assert at(b) == at(c) + at(d)
                        # Every induced restriction: absent fiber, or a smaller substitution.
                        for u in range(1 << len(s)):
                            bmask = u & fiber_mask
                            if not bmask:
                                assert faces(s, u) == faces(a, u)
                            else:
                                amask = u & old_mask
                                ka = faces(a, amask | (1 << v))
                                la = faces(a, amask & ~a[v])
                                xb = frozenset(f for f in x if f & bmask == f)
                                assert faces(s, u) == faces(a, amask) | join(la, xb)
                                assert la == frozenset(f ^ (1 << v) for f in ka if f >> v & 1)
                            restrictions += 1
                        cases += 1
    print(f'SUBSTITUTION: all labelled G of orders 1..4, H of orders 1..3, all vertices: {cases} cases')
    print(f'RESTRICTIONS: {restrictions}; empty-realization links: {link_empty}; all union/intersection/link identities pass')
    print('SPLITTING DIAGNOSTIC: rational Betti additivity passes in every nonempty-link case')


def boundary_family():
    for m in range(2, 5):
        for n in range(2, 5):
            a = graph(m, combinations(range(m), 2))
            h = graph(2 * n, [(2 * i, 2 * i + 1) for i in range(n)])
            k = faces(lex(a, h))
            # Exact equality with m disjoint boundaries of n-dimensional cross-polytopes.
            expected = {0}
            for j in range(m):
                fiber_faces = {0}
                for i in range(n):
                    fiber_faces = {f | b for f in fiber_faces for b in
                                   [0, 1 << (2 * n * j + 2 * i), 1 << (2 * n * j + 2 * i + 1)]}
                expected.update(fiber_faces)
            assert k == frozenset(expected)
            wanted = (m,) + (0,) * (n - 2) + (m,)
            assert betti(k) == wanted
            print(f'BOUNDARY: K_{m} o ({n} K_2): {len(k)-1} nonempty faces; Betti={wanted}; {m} cross-polytope boundaries S^{n-1}')


def main():
    source_fidelity()
    small_substitution_checks()
    boundary_family()
    print('LIMIT: finite combinatorial/collapse/homology checks; no universal wedge-closure theorem is certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
