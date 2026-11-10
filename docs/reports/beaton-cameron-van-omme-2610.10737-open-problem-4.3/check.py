#!/usr/bin/env python3
"""Exact finite checks for Beaton--Cameron--van Omme Open Problem 4.3.

Original implementation; standard library only. Graphs are integer adjacency
bitsets. These checks test finite instances, not the quantified settlement.
"""

from functools import lru_cache
from itertools import combinations


def bits(s):
    while s:
        b = s & -s
        yield b.bit_length() - 1
        s ^= b


def graph(n, edges):
    a = [0] * n
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def edges(g):
    return [(u, v) for u in range(len(g)) for v in bits(g[u]) if u < v]


def coloring(g, k, deleted=None):
    """Exhaustive DSATUR search, with color-permutation symmetry removed."""
    active = ((1 << len(g)) - 1) & ~(0 if deleted is None else 1 << deleted)
    c = [-1] * len(g)

    def visit(todo, used):
        if not todo:
            return tuple(c)
        forbidden = {v: sum(1 << i for i in {c[w] for w in bits(g[v]) if c[w] >= 0})
                     for v in bits(todo)}
        v = max(forbidden, key=lambda x: (forbidden[x].bit_count(), (g[x] & todo).bit_count()))
        for col in range(min(k, used + 1)):
            if not (forbidden[v] >> col) & 1:
                c[v] = col
                answer = visit(todo ^ (1 << v), max(used, col + 1))
                if answer is not None:
                    return answer
        c[v] = -1
        return None

    return visit(active, 0)


def chi(g):
    return next(k for k in range(len(g) + 1) if coloring(g, k) is not None)


def critical(g, k):
    return chi(g) == k and all(coloring(g, k - 1, v) is not None for v in range(len(g)))


def triangle_free(g):
    return all(not (g[u] & g[v]) for u, v in edges(g))


def induced_matching(g):
    """Vertex recursion: skip v, or choose one incident matching edge vw.

    Once vw is chosen, delete the closed neighborhoods of both endpoints.
    The remainder is precisely where another induced-matching edge can lie.
    """
    @lru_cache(None)
    def solve(s):
        if not s:
            return 0
        v = max(bits(s), key=lambda x: (g[x] & s).bit_count())
        neighbors = g[v] & s
        if not neighbors:
            return 0
        best = solve(s & ~(1 << v))
        for w in bits(neighbors):
            best = max(best, 1 + solve(s & ~(g[v] | g[w] | (1 << v) | (1 << w))))
        return best

    return solve((1 << len(g)) - 1)


def brute_induced_matching(g):
    """Independent subset-of-vertices check used on every small base graph."""
    n = len(g)
    for k in range(n // 2, 0, -1):
        for vertices in combinations(range(n), 2 * k):
            s = sum(1 << v for v in vertices)
            if all((g[v] & s).bit_count() == 1 for v in vertices):
                return k
    return 0


def mycielski(g):
    n = len(g)
    es = edges(g)
    return graph(2 * n + 1, es + [(u, n + v) for u, v in es]
                 + [(v, n + u) for u, v in es] + [(2 * n, n + v) for v in range(n)])


def seed(t):
    """Zhou et al. 2512.12349v1/v2, Section 4, formula preceding Theorem 4.1."""
    n = 3 * t + 10
    offsets = [1, n - 1] + [5 + 3 * j for j in range(t + 1)]
    return graph(n, ((i, (i + d) % n) for i in range(n) for d in offsets))


def paw_free(g):
    return all(sorted((g[v] & sum(1 << w for w in vs)).bit_count() for v in vs)
               != [1, 2, 2, 3] for vs in combinations(range(len(g)), 4))


def components(g):
    todo = (1 << len(g)) - 1
    while todo:
        seen = 0
        frontier = todo & -todo
        while frontier:
            seen |= frontier
            neighbors = 0
            for v in bits(frontier):
                neighbors |= g[v]
            frontier = neighbors & ~seen
        todo &= ~seen
        yield seen


def multipartite(g, s):
    # Reflexive nonadjacency is transitive exactly for complete multipartite graphs.
    return all((g[u] >> v) & 1 or (g[v] >> w) & 1 or not ((g[u] >> w) & 1)
               for u in bits(s) for v in bits(s) for w in bits(s))


def clique(g):
    return all(a.bit_count() == len(g) - 1 for a in g)


def complete_multipartite(parts):
    labels = [i for i, size in enumerate(parts) for _ in range(size)]
    return graph(len(labels), ((u, v) for u, v in combinations(range(len(labels)), 2)
                              if labels[u] != labels[v]))


def main():
    # Fidelity: the source explicitly cites the infinite (4K2,C3)-free 4-critical family.
    for t in range(3, 7):
        g = seed(t)
        assert triangle_free(g) and critical(g, 4)
        r = induced_matching(g)
        assert r <= 3
        print(f'SEED t={t} order={len(g)} edges={len(edges(g))} chi=4 critical=True nu_ind={r}')

    counts = [0, 0, 0, 0]
    for n in range(7):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            g = graph(n, (e for i, e in enumerate(pairs) if (mask >> i) & 1))
            q, r = chi(g), induced_matching(g)
            assert r == brute_induced_matching(g)
            mg = mycielski(g)
            mr = induced_matching(mg)
            assert mr <= (q + 1) * r + 1
            assert chi(mg) == q + 1
            if triangle_free(g):
                assert triangle_free(mg) and q <= 2 * r + 1
                counts[1] += 1
            if q >= 2 and critical(g, q):
                assert critical(mg, q + 1)
                counts[2] += 1
            if paw_free(g):
                for s in components(g):
                    has_triangle = any(g[u] & g[v] & s for u, v in edges(g)
                                       if (s >> u) & 1 and (s >> v) & 1)
                    if has_triangle:
                        assert multipartite(g, s)
                if critical(g, q) and not triangle_free(g):
                    assert clique(g)
                counts[3] += 1
            counts[0] += 1
        print(f'EXHAUSTIVE n={n} graphs={1 << len(pairs)} OK', flush=True)
    print(f'TOTAL graphs={counts[0]} triangle_free={counts[1]} critical_nonempty={counts[2]} paw_free={counts[3]}')

    # The criticality-preservation lemma requires q >= 2 (the route uses q >= 4).
    assert critical(graph(1, []), 1) and not critical(mycielski(graph(1, [])), 2)
    print('BOUNDARY Mycielski(K1) is not vertex-critical; lifting lemma restricted to q>=2')

    for q in range(5, 10):
        orders = []
        for t in range(1, 6):
            g = complete_multipartite([t] + [1] * (q - 1))
            assert chi(g) == q and paw_free(g) and induced_matching(g) == 1
            assert len(list(components(g))) == 1
            assert critical(g, q) == (t == 1)
            orders.append(len(g))
        print(f'HIGH_LEVEL q={q} orders={orders} chi={q} nu_ind=1 connected=True')
    r = 3
    for k in range(4, 9):
        print(f'PARAMETER k={k} r_k={r} triangle_free_cutoff={2*r+1}')
        r = (k + 1) * r + 1
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
