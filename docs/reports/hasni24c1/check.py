#!/usr/bin/env python3
"""Exact checks for the lower-bound assertion of Hasni et al., Conjecture 1.

Only Python integers, Fraction, and SymPy exact expressions are evidence.
Graph invariants are computed from adjacency, independently of the formulas.
"""
from fractions import Fraction as Q
from itertools import combinations
import sys

import sympy as sp


def graph(n, edges):
    adjacency = [set() for _ in range(n)]
    edge_set = set()
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        edge = tuple(sorted((u, v)))
        assert edge not in edge_set
        edge_set.add(edge)
        adjacency[u].add(v)
        adjacency[v].add(u)
    return adjacency, edge_set


def is_tree(adjacency, edges):
    if not adjacency:
        return False
    seen = {0}
    frontier = [0]
    while frontier:
        u = frontier.pop()
        for v in adjacency[u] - seen:
            seen.add(v)
            frontier.append(v)
    return len(seen) == len(adjacency) and len(edges) == len(adjacency) - 1


def harmonic(adjacency, edges):
    return sum((Q(2, len(adjacency[u]) + len(adjacency[v])) for u, v in edges), Q(0))


def minimum_total_domination(adjacency):
    n = len(adjacency)
    target = (1 << n) - 1
    masks = [sum(1 << v for v in ns) for ns in adjacency]
    for k in range(n + 1):
        for chosen in combinations(range(n), k):
            dominated = 0
            for u in chosen:
                dominated |= masks[u]
            if dominated == target:
                return k, chosen
    raise AssertionError('No total dominating set: graph outside source domain')


def lower_bound(n, g):
    assert 0 <= g <= n
    return (Q(n) - Q(3 * g - 1, 2)) * Q(2, n - g + 1) + Q(g - 1, 2) * (
        Q(2, n - g + 2) + Q(7, 6)
    )


def path(n):
    return graph(n, ((i, i + 1) for i in range(n - 1)))


def broom(r):
    # 0=c, 1=a, 2=b, 3=d, 4=e; 5,...,r+4 are hub leaves.
    return graph(r + 5, [(0, 1), (1, 2), (2, 3), (3, 4)] + [
        (0, i) for i in range(5, r + 5)
    ])


def source_template(r, s):
    # One hub, r leaves and s arms of length three, as in Figure 2.
    edges = [(0, i) for i in range(1, r + 1)]
    for j in range(s):
        a = r + 1 + 3 * j
        edges.extend([(0, a), (a, a + 1), (a + 1, a + 2)])
    return graph(1 + r + 3 * s, edges)


def main():
    assert sys.version_info >= (3, 12)
    print('Arithmetic: Python integers, fractions.Fraction, SymPy exact rationals')
    print('Fidelity: source Lemma 1 on its explicitly included paths P_(4t)')
    for t in range(1, 5):
        n = 4 * t
        adj, edges = path(n)
        assert is_tree(adj, edges)
        g, _ = minimum_total_domination(adj)
        h = harmonic(adj, edges)
        f = Q(11, 30) * n + Q(4, 15) * g - Q(1, 6)
        assert g == 2 * t and h == f
        print(f'P_{n}: gamma_t={g}, H={h}, f(n,g)={f}')
    adj, edges = graph(4, [(0, 1), (0, 2), (0, 3)])
    g, _ = minimum_total_domination(adj)
    assert g == 2 and harmonic(adj, edges) == Q(3, 2)
    assert Q(3, 2) < Q(11, 6)
    print('Source Theorem 1 base-case reading: H(S_4)=3/2 < f(4,2)=11/6')

    print('Independent adjacency and exhaustive total-domination check of B_r, 1 <= r <= 32')
    for r in range(1, 33):
        adj, edges = broom(r)
        n = len(adj)
        assert is_tree(adj, edges)
        assert sorted(map(len, adj)) == sorted([r + 1, 2, 2, 2] + [1] * (r + 1))
        g, witness = minimum_total_domination(adj)
        assert g == 4
        assert set(witness) == {0, 1, 2, 3}
        h = harmonic(adj, edges)
        f = lower_bound(n, g)
        assert h == Q(2 * r, r + 2) + Q(2, r + 3) + Q(5, 3)
        assert f - h == Q(1, 12) - Q(1, (r + 2) * (r + 3))
        if r >= 3:
            assert 2 * g <= n and f - h >= Q(1, 20) and h < f
        if r in (1, 2, 3, 8, 16, 32):
            print(f'B_{r}: n={n}, gamma_t={g}, H={h}, F={f}, F-H={f-h}')
    adj, edges = broom(3)
    assert harmonic(adj, edges) == Q(16, 5)
    assert lower_bound(8, 4) == Q(13, 4)

    for r in range(1, 5):
        for s in range(1, 4):
            adj, edges = source_template(r, s)
            assert is_tree(adj, edges)
            g, _ = minimum_total_domination(adj)
            assert g == 2 * s + 1
            assert harmonic(adj, edges) == lower_bound(len(adj), g)
    print('Figure 2 template: exact equality and odd gamma_t checked for 1<=r<=4, 1<=s<=3')

    r = sp.Symbol('r', integer=True, positive=True)
    n, g = r + 5, sp.Integer(4)
    f = (n - sp.Rational(1, 2) * (3 * g - 1)) * 2 / (n - g + 1) + (
        (g - 1) / 2 * (2 / (n - g + 2) + sp.Rational(7, 6))
    )
    h = 2 * r / (r + 2) + 2 / (r + 3) + sp.Rational(5, 3)
    gap = sp.Rational(1, 12) - 1 / ((r + 2) * (r + 3))
    assert not any(x.has(sp.Float) for x in (f, h, gap))
    assert sp.cancel(f - h - gap) == 0
    assert sp.cancel(gap - sp.Rational(1, 20) - (
        (r - 3) * (r + 8) / (30 * (r + 2) * (r + 3))
    )) == 0
    k = sp.Symbol('k', nonnegative=True, integer=True)
    num, den = sp.fraction(sp.factor((gap - sp.Rational(1, 20)).subs(r, k + 3)))
    assert all(c >= 0 for c in sp.Poly(num, k).all_coeffs())
    assert all(c >= 0 for c in sp.Poly(den, k).all_coeffs()) and den.subs(k, 0) > 0
    assert sp.limit(gap, r, sp.oo) == sp.Rational(1, 12)
    print('Symbolic identity: F-H = 1/12 - 1/((r+2)(r+3))')
    print('Nonnegative-coefficient certificate at r=k+3: F-H-1/20 =', sp.factor(num / den))
    print('Exact limit of F-H as r tends to infinity: 1/12')
    print('Scope: finite graph checks plus symbolic arithmetic; no Lean settlement/probe')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
