#!/usr/bin/env python3
"""Exact finite checks for Buckley--Seco, arXiv:2610.40089v1.

Finite matrices test identities, never the spectrum of an infinite operator.
Only integers and fractions.Fraction are used as mathematical evidence.
"""
from fractions import Fraction as Q
from itertools import product
from math import factorial


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum((x * y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def inverse(a):
    """Gauss--Jordan elimination, independent of the ancestor formula."""
    n = len(a)
    b = [row[:] + unit for row, unit in zip(a, eye(n))]
    for j in range(n):
        k = next(k for k in range(j, n) if b[k][j])
        b[j], b[k] = b[k], b[j]
        pivot = b[j][j]
        b[j] = [x / pivot for x in b[j]]
        for k in range(n):
            if k != j:
                c = b[k][j]
                b[k] = [x - c * y for x, y in zip(b[k], b[j])]
    return [row[n:] for row in b]


def ancestor(parent, v, m):
    for _ in range(m):
        v = parent[v]
        if v is None:
            return None
    return v


def distance_below(parent, w, v):
    for m in range(len(parent)):
        if w == v:
            return m
        w = parent[w]
        if w is None:
            return None
    raise AssertionError("cycle in a test tree")


def check_all_small():
    trees = fidelity = columns = 0
    for n in range(1, 8):
        # Every recursively labelled rooted tree: parent(j) < j, root 0.
        for choices in product(*(range(j) for j in range(1, n))):
            parent = [None, *choices]
            trees += 1
            s = [[Q(parent[j] == i) for j in range(n)] for i in range(n)]
            t = transpose(s)  # S*: coefficients pull back along the parent.
            sm, tm = eye(n), eye(n)
            for m in range(n + 1):
                counts = [sum(ancestor(parent, w, m) == v for w in range(n))
                          for v in range(n)]
                gram = mul(sm, tm)
                assert gram == [[Q(counts[i] if i == j else 0)
                                 for j in range(n)] for i in range(n)]
                # A diagonal nonnegative Gram matrix has this exact top eigenvalue.
                assert max(gram[i][i] for i in range(n)) == max(counts)
                fidelity += 1
                sm, tm = mul(sm, s), mul(tm, t)
            for z in (Q(5, 4), Q(-3, 2)):
                a = [[t[i][j] - z * (i == j) for j in range(n)] for i in range(n)]
                r = inverse(a)
                assert mul(a, r) == eye(n)
                for v in range(n):
                    norm2 = sum((r[w][v] ** 2 for w in range(n)), Q(0))
                    for w in range(n):
                        m = distance_below(parent, w, v)
                        expected = Q(0) if m is None else -z ** (-m - 1)
                        assert r[w][v] == expected
                    for m in range(n + 1):
                        count = sum(ancestor(parent, w, m) == v for w in range(n))
                        assert count * z ** (-2 * m - 2) <= norm2
                    columns += 1
    assert trees == sum(factorial(n - 1) for n in range(1, 8)) == 874
    print(f"SMALL_TREES={trees}; GRAM_IDENTITIES={fidelity}; INVERSE_COLUMNS={columns}")


def family_children(v):
    tag, label = v
    if tag == "a":
        return [("o", 0)] if label == 1 else [("a", label - 1)]
    if tag == "o":
        return [("u", 0), ("b", "")]
    if tag == "u":
        return [("u", label + 1)]
    return [("b", label + "0"), ("b", label + "1")]


def descendants_at(v, m):
    generation = [v]
    for _ in range(m):
        generation = [w for u in generation for w in family_children(u)]
    return generation


def family_generation_size(v, m):
    tag, label = v
    if tag == "u":
        return 1
    if tag == "b":
        return 2 ** m
    if tag == "o":
        return 1 if m == 0 else 1 + 2 ** (m - 1)
    return 1 if m <= label else 1 + 2 ** (m - label - 1)


def check_family():
    # Rooted example omits all a_j. Rootless example includes every a_j, j >= 1.
    seeds = [("o", 0), ("u", 0), ("b", ""), ("b", "10")]
    seeds += [("a", j) for j in range(1, 7)]
    tests = 0
    for v in seeds:
        for m in range(9):
            gen = descendants_at(v, m)
            assert len(gen) == len(set(gen)) == family_generation_size(v, m)
            assert 1 <= len(gen) <= 2 ** m
            tests += 1
    # The unary and binary vertices attain the global extrema at each depth.
    for m in range(81):
        assert family_generation_size(("u", 0), m) == 1
        assert family_generation_size(("b", ""), m) == 2 ** m
    # Prefix-free branch vertices give the source's pairwise-disjoint family A4.
    prefix_free = ["1" * n + "0" for n in range(30)]
    assert all(not a.startswith(b) and not b.startswith(a)
               for i, a in enumerate(prefix_free) for b in prefix_free[i + 1:])
    # Truncated rootless ancestor window: solve the triangular system directly.
    top, depth = ("a", 6), 13
    levels = [descendants_at(top, m) for m in range(depth + 1)]
    parent = {w: u for level in levels[:-1] for u in level for w in family_children(u)}
    target, z = ("b", ""), Q(5, 4)
    x = {}
    for level in levels:
        for w in level:
            forcing = Q(w == target)
            x[w] = (x.get(parent.get(w), Q(0)) - forcing) / z
            assert x.get(parent.get(w), Q(0)) - z * x[w] == forcing
    assert all(x[("a", j)] == 0 for j in range(1, 7))
    for m in range(7):
        gen = descendants_at(target, m)
        assert all(x[w] == -z ** (-m - 1) for w in gen)
        assert sum((x[w] ** 2 for w in gen), Q(0)) == 2 ** m * z ** (-2 * m - 2)
    # Exact geometric growth for a radius strictly inside (1, sqrt(2)).
    assert 1 < z and z ** 2 < 2
    lower = [Q(2 ** m) * z ** (-2 * m - 2) for m in range(81)]
    assert all(lower[m + 1] == Q(32, 25) * lower[m] for m in range(80))
    assert lower[80] > 100000000
    print(f"FAMILY_GENERATION_ENUMERATIONS={tests}; EXACT_EXTREMA_DEPTHS=0..80")
    print("FAMILY_RADII_SQUARED: R_minus^2=1; R_plus^2=2 (from exact growth profiles)")
    print("ROOTLESS_WINDOW: ancestor coordinates vanish; descendant formula and residual hold")
    print("INTERIOR_RADIUS z=5/4: generation norm-square ratio=32/25; depth-80 lower bound>10^8")
    print("LIMIT: finite checks support identities; no infinite spectrum is inferred from truncations")


if __name__ == "__main__":
    check_all_small()
    check_family()
    print("ALL_CHECKS_OK")
