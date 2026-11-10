#!/usr/bin/env python3
"""Exact independent checks of Lu--Li Question 8.6's final necessity clause.

Standard library only. All linear algebra is over fractions.Fraction.
The off-degree identity minors exclude EVERY real eigenvalue simultaneously.
"""
from collections import deque
from fractions import Fraction as Q
from itertools import combinations
from math import prod


def seed(extra=False):
    edges = [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5), (4, 5),
             (3, 6), (3, 7), (6, 7), (3, 8)]
    if extra:
        edges.append((0, 3))
    return graph(9, edges)


def graph(n, edges):
    g = [set() for _ in range(n)]
    for u, v in edges:
        assert u != v and 0 <= u < n and 0 <= v < n
        g[u].add(v)
        g[v].add(u)
    return g


def corona(g):
    n = len(g)
    return graph(2 * n, [(u, v) for u in range(n) for v in g[u] if u < v]
                 + [(u, n + u) for u in range(n)])


def apply_a(g, v):
    return tuple(sum((v[j] for j in g[i]), Q(0)) for i in range(len(g)))


def apply_d(g, v):
    return tuple(len(g[i]) * v[i] for i in range(len(g)))


class Span:
    def __init__(self, n):
        self.n, self.rows = n, {}

    def add(self, v):
        w = list(map(Q, v))
        assert len(w) == self.n
        for pivot, row in sorted(self.rows.items()):
            c = w[pivot]
            if c:
                w = [a - c * b for a, b in zip(w, row)]
        if not any(w):
            return None
        pivot = next(i for i, x in enumerate(w) if x)
        c = w[pivot]
        w = tuple(x / c for x in w)
        self.rows[pivot] = w
        return w

    def contains(self, v):
        w = list(map(Q, v))
        for pivot, row in sorted(self.rows.items()):
            c = w[pivot]
            w = [a - c * b for a, b in zip(w, row)]
        return not any(w)


def cyclic_space(g):
    s = Span(len(g))
    todo = deque([s.add([1] * len(g))])
    while todo:
        v = todo.popleft()
        for w in (apply_a(g, v), apply_d(g, v)):
            new = s.add(w)
            if new is not None:
                todo.append(new)
    # Independent closure verification: basis is invariant and contains 1.
    assert s.contains([1] * len(g))
    for v in s.rows.values():
        assert s.contains(apply_a(g, v)) and s.contains(apply_d(g, v))
    return s


def stable_cells(g):
    colors = [0] * len(g)
    while True:
        signatures = [(colors[v], tuple(sorted(colors[u] for u in g[v])))
                      for v in range(len(g))]
        table = {s: i for i, s in enumerate(sorted(set(signatures)))}
        new = [table[s] for s in signatures]
        old_cells = {frozenset(i for i, c in enumerate(colors) if c == x)
                     for x in set(colors)}
        cells = {frozenset(i for i, c in enumerate(new) if c == x)
                 for x in set(new)}
        if cells == old_cells:
            # Every final cell has equal neighbour counts into every cell.
            for c in cells:
                for d in cells:
                    assert len({len(g[u] & d) for u in c}) == 1
            return cells
        colors = new


def bfs(g, source, allowed):
    dist = {source: 0}
    todo = deque([source])
    while todo:
        v = todo.popleft()
        for w in g[v] & allowed:
            if w not in dist:
                dist[w] = dist[v] + 1
                todo.append(w)
    return dist


def exhaustive_dh(g):
    n = len(g)
    whole = set(range(n))
    ambient = [bfs(g, v, whole) for v in range(n)]
    assert len(ambient[0]) == n
    checked = 0
    for mask in range(1, 1 << n):
        s = {v for v in range(n) if (mask >> v) & 1}
        first = bfs(g, min(s), s)
        if len(first) != len(s):
            continue
        checked += 1
        for u in s:
            assert bfs(g, u, s) == {v: ambient[u][v] for v in s}
    return checked


def dh_reduction(g):
    # Source's pendant/twin characterization, with each local step verified.
    active, trace = set(range(len(g))), []
    while len(active) > 1:
        for u in sorted(active):
            nbr = g[u] & active
            if len(nbr) == 1:
                trace.append(('pendant', u, next(iter(nbr))))
                active.remove(u)
                break
        else:
            for u, v in combinations(sorted(active), 2):
                if (g[u] & active) - {v} == (g[v] & active) - {u}:
                    trace.append(('true-twin' if v in g[u] else 'false-twin', u, v))
                    active.remove(u)
                    break
            else:
                raise AssertionError('No pendant/twin reduction')
    return trace


def no_degree_eigenvectors(g, parent_size):
    n = parent_size
    assert len(g) == 2 * n
    assert all(g[n + u] == {u} and len(g[u]) > 1 for u in range(n))
    for d in sorted({len(x) for x in g}):
        cols = sorted(u for u in range(2 * n) if len(g[u]) == d)
        rows = [u - n if d == 1 else u + n for u in cols]
        assert all(len(g[r]) != d for r in rows)
        minor = [[int(c in g[r]) for c in cols] for r in rows]
        assert minor == [[int(i == j) for j in range(len(cols))]
                         for i in range(len(cols))]
        print(f'  degree={d}: columns={cols}, outside_rows={rows}, '
              f'minor=I_{len(cols)}, determinant=1; all real lambda excluded')


def main():
    f = (0, -2, -2, 0, 1, 1, 1, 1, 0)
    source_cells = {frozenset([4, 5]), frozenset([6, 7])}
    source_cells |= {frozenset([v]) for v in [0, 1, 2, 3, 8]}
    for extra in [False, True]:
        g = seed(extra)
        s, cells = cyclic_space(g), stable_cells(g)
        assert len(s.rows) == 6 and cells == source_cells
        assert apply_a(g, f) == f and apply_d(g, f) == tuple(2 * x for x in f)
        assert sum(f) == 0 and not s.contains(f)
        assert all(sum(v[i] * f[i] for i in range(9)) == 0 for v in s.rows.values())
        walk_span, v = Span(9), tuple([1] * 9)
        for _ in range(6):
            walk_span.add(v)
            v = apply_a(g, v)
        assert len(walk_span.rows) == 6
        assert all(len({f[i] for i in c}) == 1 for c in cells)
        print(f'FIDELITY G{2 if extra else 1}: n=9, dim M=6, dim W=7; '
              'Af=f, Df=2f, sum f=0, f in W and perpendicular to M; walk rank=6')
        print('  DH reduction:', dh_reduction(g))
        print('  connected induced subsets checked:', exhaustive_dh(g))
        for t in [1, 2]:
            old_g, old_s, old_cells = g, s, cells
            g = corona(old_g)
            s, cells = cyclic_space(g), stable_cells(g)
            n = len(old_g)
            expected_cells = old_cells | {frozenset(n + u for u in c) for c in old_cells}
            assert cells == expected_cells
            doubled = Span(2 * n)
            for v in old_s.rows.values():
                doubled.add(v + (Q(0),) * n)
                doubled.add((Q(0),) * n + v)
            assert len(s.rows) == len(doubled.rows) == 6 * 2**t
            assert len(cells) == 7 * 2**t
            assert all(s.contains(v) for v in doubled.rows.values())
            assert all(doubled.contains(v) for v in s.rows.values())
            degrees = sorted({len(x) for x in old_g})
            for u in range(2 * n):
                q = prod(Q(len(g[u]) - d - 1, -d) for d in degrees)
                assert q == int(u >= n)
            print(f'COUNTEREXAMPLE G{2 if extra else 1} corona^{t}: '
                  f'n={2*n}, edges={sum(map(len,g))//2}, '
                  f'dim M={len(s.rows)}, dim W={len(cells)}; '
                  'M and W equal the two lifted copies; leaf polynomial verified')
            print('  DH reduction:', dh_reduction(g))
            if t == 1:
                print('  connected induced subsets checked:', exhaustive_dh(g))
            no_degree_eigenvectors(g, n)
            if not extra and t == 1:
                u, v = f + (0,) * 9, (0,) * 9 + f
                assert apply_a(g, u) == tuple(a+b for a,b in zip(u,v))
                assert apply_a(g, v) == u
                assert apply_d(g, u) == tuple(3*x for x in u)
                assert apply_d(g, v) == v
                assert all(sum(w[i]*z[i] for i in range(18)) == 0
                           for w in s.rows.values() for z in [u,v])
                print('  hidden plane: A=[[1,1],[1,0]], D=diag(3,1), orthogonal to M')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
