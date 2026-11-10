#!/usr/bin/env python3
"""Exact preregistration checks for arXiv:2610.14681v1, Problem 4.

Only Python integers and finite sets are used. This is independent experimental
support for the proposed argument, not a Lean proof or an unbounded enumeration.
"""
from collections import deque
from itertools import combinations, product


def cycle(n):
    return {i: {(i - 1) % n, (i + 1) % n} for i in range(n)}


def box(h, k):
    return {(a, b): {(x, b) for x in h[a]} | {(a, y) for y in k[b]}
            for a in h for b in k}


def source_labels(n):
    # The source theorem labelled theo:cycles, including its even and odd displayed formulae.
    if n == 3:
        return [1, 1, 1], 0, 0
    k = n // 2
    if n % 2 == 0:
        p = [1 if i == 0 else 4 * i if i < k else 4 * k + 1 if i == k
             else 4 * (2 * k - i) + 2 for i in range(n)]
    else:
        p = [1 if i == 0 else 4 * i if i <= k else 4 * k + 3 if i == k + 1
             else 4 * (2 * k + 1 - i) + 2 for i in range(n)]
    return p, 3, 5


def distances(g, start):
    d, q = {start: 0}, deque([start])
    while q:
        v = q.popleft()
        for w in g[v] - d.keys():
            d[w] = d[v] + 1
            q.append(w)
    return d


def fidelity():
    pair_checks = 0
    for n in range(3, 33):
        p, a, b = source_labels(n)
        g = cycle(n)
        assert min(p) > 0
        # Construct the actual unit-edge caterpillar, then measure by BFS.
        r = max(p)
        t = {i: set() for i in range(r + n)}
        for i in range(r - 1):
            t[i].add(i + 1)
            t[i + 1].add(i)
        for i, label in enumerate(p):
            t[r + i].add(label - 1)
            t[label - 1].add(r + i)
        assert len(distances(t, 0)) == len(t)
        assert sum(map(len, t.values())) == 2 * (len(t) - 1)
        leaves = {v for v in t if len(t[v]) == 1}
        assert leaves == set(range(r, r + n))
        core = set(t) - leaves
        assert core == set(range(r))
        assert all(t[v] & core == ({v - 1, v + 1} & core) for v in core)
        for i in g:
            d = distances(t, r + i)
            for j in g:
                if i == j:
                    continue
                delta = abs(p[i] - p[j])
                assert (j in g[i]) == (a <= delta <= b)
                assert d[r + j] == delta + 2
                assert (j in g[i]) == (a + 2 <= d[r + j] <= b + 2)
                pair_checks += 1
    print(f'FIDELITY: C_n for 3 <= n <= 32, labels and literal caterpillar distances: {pair_checks} ordered pairs OK')


def scalar_check():
    cases = 0
    # Translate the global minimum to 0. x and z are the extremal neighbor labels.
    for a in range(7):
        for b in range(a, 9):
            for x in range(a, b + 1):
                for z in range(x, b + 1):
                    for y in range(x, z + 1):
                        for w in range(z + b + 1):
                            if a <= abs(w - y) <= b:
                                assert a <= abs(w - x) <= b or a <= abs(w - z) <= b
                                cases += 1
    print(f'SCALAR: 0 <= a <= 6, a <= b <= 8, all in-range x <= y <= z and w: {cases} premises OK')


def obstruction_check(g):
    # Every possible pair of extremal neighbors is defeated. This is the
    # combinatorial obstruction, independently of any proposed label values.
    count = 0
    for v in g:
        for x, z in combinations(sorted(g[v]), 2):
            witnesses = [(y, w) for y in g[v] - {x, z} for w in g[y] - {x, z}
                         if w not in g[x] and w not in g[z]]
            assert witnesses, (v, x, z)
            for y, w in witnesses:
                assert y != x and y != z and w != x and w != z
                assert y in g[v] and w in g[y]
            count += 1
    return count


def named_construction_check(h, k):
    g = box(h, k)
    count = 0
    for v in g:
        h0, k0 = v
        for x, z in combinations(sorted(g[v]), 2):
            if x[1] == k0 and z[1] == k0:
                k1 = min(k[k0])
                k2 = min(k[k1] - {k0})
                y, w = (h0, k1), (h0, k2)
            elif x[0] == h0 and z[0] == h0:
                h1 = min(h[h0])
                h2 = min(h[h1] - {h0})
                y, w = (h1, k0), (h2, k0)
            else:
                hx = x if x[1] == k0 else z
                kz = z if x[1] == k0 else x
                h2 = min(h[h0] - {hx[0]})
                k2 = min(k[k0] - {kz[1]})
                y, w = (h2, k0), (h2, k2)
            assert y in g[v] - {x, z}
            assert w in g[y] - {x, z}
            assert w not in g[x] and w not in g[z]
            count += 1
    assert count == obstruction_check(g)
    return count


def settlement_check():
    count = sum(named_construction_check(cycle(m), cycle(n))
                for m, n in product(range(3, 13), repeat=2))
    print(f'TOROIDAL: all 100 ordered (m,n) with 3 <= m,n <= 12: {count} vertex/extreme-pair obstructions OK')
    complete = lambda n: {v: set(range(n)) - {v} for v in range(n)}
    bipartite23 = {v: set(range(2, 5)) if v < 2 else {0, 1} for v in range(5)}
    disconnected = {**cycle(3), **{v + 3: {w + 3 for w in ns} for v, ns in cycle(3).items()}}
    factors = [cycle(3), cycle(4), complete(4), complete(5), bipartite23, disconnected]
    extra = sum(named_construction_check(h, k) for h, k in product(factors, repeat=2))
    print(f'GENERAL PRODUCTS: 36 ordered pairs from C3,C4,K4,K5,K2,3,2C3: {extra} obstructions OK')
    # Negative control: the obstruction must not reject a cycle, which the source represents.
    for n in range(3, 13):
        g = cycle(n)
        good = []
        for v in g:
            x, z = sorted(g[v])
            good.append(not any(w not in g[x] and w not in g[z]
                                for y in g[v] - {x, z} for w in g[y] - {x, z}))
        assert all(good)
    print('NEGATIVE CONTROL: represented cycles do not trigger the obstruction')


def main():
    fidelity()
    scalar_check()
    settlement_check()
    print('SCOPE: exact finite checks supporting a general argument; no unbounded theorem certified here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
