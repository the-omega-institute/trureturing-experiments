#!/usr/bin/env python3
"""Exact finite checks for Sason--Krupnik--Hamud--Berman Question 4.

Graphs are built from edge sets, independently of the proposed block transport.
Requires Python >= 3.12 and SymPy. All evidence is integer/rational arithmetic.
"""

from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
import platform

import sympy as sp


def graph(n, edges):
    edges = tuple(sorted({tuple(sorted(e)) for e in edges}))
    assert all(0 <= u < v < n for u, v in edges)
    return n, edges


def cycle(n):
    return graph(n, ((i, (i + 1) % n) for i in range(n)))


def complete(n):
    return graph(n, combinations(range(n), 2))


def path(n):
    return graph(n, ((i, i + 1) for i in range(n - 1)))


@lru_cache(None)
def adjacency(g):
    n, edges = g
    a = sp.zeros(n)
    for u, v in edges:
        a[u, v] = a[v, u] = 1
    return sp.ImmutableMatrix(a)


def degrees(g):
    return tuple(sum(adjacency(g).row(i)) for i in range(g[0]))


@lru_cache(None)
def coefficients(g, kind="Q"):
    a = adjacency(g)
    d = sp.diag(*degrees(g)) if g[0] else sp.zeros(0)
    if kind == "A":
        m = a
    elif kind == "L":
        m = d - a
    elif kind == "Q":
        m = d + a
    elif kind == "N":
        # No isolates in these fidelity cases. This is similar to D^-1/2 L D^-1/2.
        assert all(t > 0 for t in degrees(g))
        m = sp.eye(g[0]) - sp.diag(*(sp.Rational(1, t) for t in degrees(g))) * a
    else:
        raise ValueError(kind)
    result = tuple(m.charpoly().all_coeffs())
    assert all(c.is_Rational for c in result)
    return result


def expansion(g, retain):
    n, edges = g
    result = list(edges) if retain else []
    for k, (u, v) in enumerate(edges):
        result.extend(((u, n + k), (v, n + k)))
    return graph(n + len(edges), result)


def partial_join(g, h, left_inserted, right_inserted, retain=True):
    """Literal disjoint expansions, then edges between the selected classes."""
    s, b = expansion(g, False), expansion(h, retain)
    offset = s[0]
    edges = list(s[1]) + [(u + offset, v + offset) for u, v in b[1]]
    left = range(g[0], s[0]) if left_inserted else range(g[0])
    right = range(h[0], b[0]) if right_inserted else range(h[0])
    edges.extend((u, v + offset) for u in left for v in right)
    result = graph(s[0] + b[0], edges)
    base = 2 * len(g[1]) + (3 if retain else 2) * len(h[1])
    assert len(result[1]) == base + len(left) * len(right)
    return result


def regular_graphs_through_four():
    result = []
    for n in range(5):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            g = graph(n, (e for k, e in enumerate(possible) if mask & (1 << k)))
            if len(set(degrees(g))) <= 1:
                result.append(g)
    return result


def report(label, **fields):
    print(json.dumps({"check": label, **fields}, sort_keys=True), flush=True)


def main():
    report("environment", python=platform.python_version(), sympy=sp.__version__,
           arithmetic="exact SymPy integers and rationals")
    x = sp.Symbol("lambda")
    c4, p3 = cycle(4), path(3)
    sc4 = expansion(c4, False)
    expected = sp.Poly(x**2 * (x**2 - 4) * (x**2 - 2)**2, x)
    assert coefficients(sc4, "A") == tuple(expected.all_coeffs())
    report("source Definition 42: S(C4) is C8", vertices=sc4[0], edges=len(sc4[1]),
           adjacency_charpoly=str(expected.as_expr()))
    for i, j in product((False, True), repeat=2):
        joined = partial_join(c4, p3, i, j)
        assert joined[0] == 13
        expected_edges = 14 + (4 if i else 4) * (2 if j else 3)
        assert len(joined[1]) == expected_edges
        report("source Definitions 44--47: C4,P3", inserted=[i, j],
               vertices=joined[0], edges=len(joined[1]))

    small = regular_graphs_through_four()
    groups = defaultdict(list)
    for g in small:
        groups[coefficients(g, "A")].append(g)
    matrices = 0
    for retain, i, j in product((False, True), repeat=3):
        for left_group, right_group in product(groups.values(), repeat=2):
            seen = {coefficients(partial_join(g, h, i, j, retain))
                    for g, h in product(left_group, right_group)}
            assert len(seen) == 1
            matrices += len(left_group) * len(right_group)
    report("all labeled regular inputs of orders 0--4", graphs=len(small),
           adjacency_spectral_classes=len(groups), q_matrices=matrices,
           operations=4, readings=2, conclusion="one Q polynomial per pair of input spectral classes")

    vertices = list(product(range(4), repeat=2))
    rook = graph(16, ((u, v) for u, v in combinations(range(16), 2)
                     if any(vertices[u][k] == vertices[v][k] for k in (0, 1))))
    differences = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    shrikhande = graph(16, ((u, v) for u, v in combinations(range(16), 2)
                           if tuple((vertices[v][k] - vertices[u][k]) % 4
                                    for k in (0, 1)) in differences))
    assert set(degrees(rook)) == set(degrees(shrikhande)) == {6}
    assert coefficients(rook, "A") == coefficients(shrikhande, "A")
    edge_sets = [set(g[1]) for g in (rook, shrikhande)]
    k4_counts = [sum(all(e in es for e in combinations(s, 2))
                     for s in combinations(range(16), 4)) for es in edge_sets]
    assert k4_counts == [8, 0]
    report("nonisomorphic regular cospectral seeds", vertices=16, degree=6,
           edges=48, k4_counts=k4_counts,
           adjacency_coefficients=list(map(str, coefficients(rook, "A"))))

    k3 = complete(3)
    for i, j in product((False, True), repeat=2):
        left = partial_join(rook, k3, i, j)
        right = partial_join(shrikhande, k3, i, j)
        assert left[0] == right[0] == 70
        for kind in ("A", "L", "N", "Q"):
            a, b = coefficients(left, kind), coefficients(right, kind)
            assert a == b
            report("Theorem 41 fidelity" if kind != "Q" else "Question 4 expectation",
                   inserted=[i, j], reading="retained original edges", matrix=kind,
                   vertices=70, coefficients=list(map(str, a)))
        a = partial_join(rook, k3, i, j, retain=False)
        b = partial_join(shrikhande, k3, i, j, retain=False)
        assert coefficients(a) == coefficients(b)
        report("Question 4 alternate reading", inserted=[i, j], vertices=70,
               coefficients=list(map(str, coefficients(a))))
        # Replace the second input pair as well, independently of the first position.
        a = partial_join(k3, rook, i, j)
        b = partial_join(k3, shrikhande, i, j)
        assert coefficients(a) == coefficients(b)
        report("Question 4 second input varies", inserted=[i, j], vertices=70,
               coefficient_sha256=sha256(str(coefficients(a)).encode()).hexdigest())

    for retain in (False, True):
        both = [partial_join(g, h, False, False, retain)
                for g, h in product((rook, shrikhande), repeat=2)]
        assert len({coefficients(g) for g in both}) == 1
        report("both nonisomorphic input pairs vary", inserted=[False, False],
               retain=retain, vertices=128, combinations=4,
               coefficient_sha256=sha256(str(coefficients(both[0])).encode()).hexdigest())
    print("ALL_CHECKS_OK", flush=True)


if __name__ == "__main__":
    main()
