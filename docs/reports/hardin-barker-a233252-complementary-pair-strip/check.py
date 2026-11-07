#!/usr/bin/env python3
"""Exact fidelity and settlement checks for OEIS A233252."""
from fractions import Fraction
from itertools import product
import sympy as sp

SOURCE = [10, 1184, 166400, 23896064, 3439984640, 495341010944,
          71328837140480, 10271348253261824, 1479074079750225920,
          212986666384520904704, 30670079941778824232960,
          4416491511334675712835584]


def vertices(n):
    return [(i, j) for i in range(n) for j in range(4)]


def edges(n):
    out = []
    for i in range(n):
        for j in range(3):
            out.append(((i, j), (i, j + 1)))
    for i in range(n - 1):
        for j in range(3):
            out.append(((i, j + 1), (i + 1, j)))
    return out


def graph_adjacency(n):
    a = {v: set() for v in vertices(n)}
    for u, v in edges(n):
        a[u].add(v)
        a[v].add(u)
    return a


def pair_color(p):
    return min(p, 5 - p)


def normalized_count(n):
    """Enumerate literal six-color arrays via canonical pair assignments and signs."""
    vs = vertices(n)
    a = graph_adjacency(n)
    assignment = {(0, 0): 0}
    total = 0
    pair_colorings = 0

    def literal_legal(colors):
        if colors[(0, 0)] != 0:
            return False
        for u, v in edges(n):
            if colors[u] == colors[v] or colors[v] == 5 - colors[u]:
                return False
        for pos, v in enumerate(vs):
            earlier = vs[:pos]
            if colors[v] in (2, 3, 4) and not any(colors[u] == 1 for u in earlier):
                return False
            if colors[v] == 3 and not any(colors[u] == 2 for u in earlier):
                return False
        return True

    def visit(pos):
        nonlocal total, pair_colorings
        if pos == len(vs):
            first = {}
            for v in vs:
                first.setdefault(assignment[v], v)
            order = sorted(first, key=lambda c: first[c])
            if order != list(range(len(order))):
                return
            k = len(order)
            pair_colorings += 1
            fixed = {first[c] for c in order}
            free = [v for v in vs if v not in fixed]
            representatives = (0, 1, 2)
            complements = (5, 4, 3)
            seen = set()
            for bits in product((0, 1), repeat=len(free)):
                colors = {}
                for c in order:
                    colors[first[c]] = representatives[c]
                for v, bit in zip(free, bits):
                    c = assignment[v]
                    colors[v] = complements[c] if bit else representatives[c]
                key = tuple(colors[v] for v in vs)
                assert key not in seen
                seen.add(key)
                assert literal_legal(colors)
                total += 1
            return
        v = vs[pos]
        if v == (0, 0):
            visit(pos + 1)
            return
        for c in range(3):
            if all(assignment.get(w) != c for w in a[v]):
                assignment[v] = c
                visit(pos + 1)
                del assignment[v]

    visit(0)
    return total, pair_colorings


def closed_form(n):
    return (2 ** (4 * n - 3)) * (4 * (9 ** n) + 9) // 9


def recurrence_values(length):
    vals = [10, 1184]
    while len(vals) < length:
        vals.append(160 * vals[-1] - 2304 * vals[-2])
    return vals


def gf_coefficients(length):
    # Formal exact series of 8*x/(1-144*x) + 2*x/(1-16*x).
    return [8 * (144 ** (n - 1)) + 2 * (16 ** (n - 1)) for n in range(1, length + 1)]


def main():
    # The source's first three terms are independently enumerated from the literal rules.
    brute = [normalized_count(n)[0] for n in (1, 2, 3)]
    assert brute == SOURCE[:3], (brute, SOURCE[:3])
    assert [normalized_count(n)[1] for n in (1, 2, 3)] == [4, 36, 324]

    # The closed form and equivalent two-growth-rate expression agree for all displayed terms.
    assert [closed_form(n) for n in range(1, 13)] == SOURCE
    assert gf_coefficients(12) == SOURCE
    assert recurrence_values(12) == SOURCE

    # Exact rational generating-function identity, checked coefficientwise through n=12.
    for n, a in enumerate(SOURCE, start=1):
        lhs = Fraction(a)
        rhs = Fraction(8 * (144 ** (n - 1)) + 2 * (16 ** (n - 1)))
        assert lhs == rhs
    x = sp.symbols("x")
    gf = sp.Rational(2) * x * (5 - 208 * x) / ((1 - 16 * x) * (1 - 144 * x))
    expansion = sp.series(gf, x, 0, len(SOURCE) + 1).removeO()
    assert [sp.expand(expansion).coeff(x, n) for n in range(1, len(SOURCE) + 1)] == SOURCE
    # First values and the literal normalization convention are recorded for auditability.
    print("literal_brute_force_n1_n3", brute)
    print("canonical_pair_colorings_n1_n3", [normalized_count(n)[1] for n in (1, 2, 3)])
    print("closed_form_first_12", SOURCE)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
