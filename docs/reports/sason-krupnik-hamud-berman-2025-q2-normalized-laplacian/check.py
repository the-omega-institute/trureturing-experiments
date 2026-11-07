#!/usr/bin/env python3
"""Exact diagnostics for SKHB Question 2 / arXiv:2412.20775v6 Question 5.19.

No finite check proves the universal claim. Matrices are integers and Fractions;
no floating point is used. SymPy supplies integer characteristic polynomials.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import lcm
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def graph(n, edges):
    a = [[0] * n for _ in range(n)]
    for i, j in edges:
        assert i != j
        a[i][j] = a[j][i] = 1
    return tuple(map(tuple, a))


def degree(a):
    return tuple(map(sum, a))


@lru_cache(None)
def integer_charpoly(a):
    n = len(a)
    if not n:
        return (1,)
    return tuple(map(int, DomainMatrix.from_list_sympy(n, n, list(map(list, a))).charpoly()))


@lru_cache(None)
def charpoly(a, kind="A"):
    n = len(a)
    d = degree(a)
    if kind == "A":
        return integer_charpoly(a)
    if kind in ("L", "Q"):
        sign = -1 if kind == "L" else 1
        return integer_charpoly(tuple(tuple((d[i] if i == j else 0) + sign * a[i][j]
                                            for j in range(n)) for i in range(n)))
    assert kind == "NL"
    # The isolated-vertex row is ZERO. On positive-degree vertices this is
    # D^(-1)(D-A), similar to D^(-1/2)(D-A)D^(-1/2). Isolated vertices split off.
    scale = lcm(*(x for x in d if x)) if any(d) else 1
    integer_matrix = tuple(tuple((scale if i == j else 0) - scale // d[i] * a[i][j]
                                 if d[i] else 0 for j in range(n)) for i in range(n))
    return tuple(Fraction(c, scale ** i) for i, c in enumerate(integer_charpoly(integer_matrix)))


@lru_cache(None)
def corona(a, b, kind):
    n, h = len(a), len(b)
    original_edges = [(i, j) for i in range(n) for j in range(i + 1, n) if a[i][j]]
    indices = original_edges if kind == "edge" else list(range(n))
    size = 2 * n + len(indices) * h
    edges = [(i, n + j) for i in range(n) for j in range(n) if a[i][j]]
    for k, index in enumerate(indices):
        offset = 2 * n + k * h
        edges += [(offset + i, offset + j) for i in range(h) for j in range(i + 1, h) if b[i][j]]
        if kind == "corona":
            endpoints = [index]
        elif kind == "neighborhood":
            endpoints = [n + j for j in range(n) if a[index][j]]
        else:
            assert kind == "edge"
            endpoints = list(index)  # ORIGINAL endpoints, even though their edge is absent in Du(a).
        edges += [(v, offset + j) for v in endpoints for j in range(h)]
    result = graph(size, edges)
    r = degree(a)[0] if n else 0
    s = degree(b)[0] if h else 0
    if len(set(degree(a))) <= 1 and len(set(degree(b))) <= 1:
        expected = {"corona": (r + h, r, s + 1),
                    "neighborhood": (r, r * (1 + h), s + r),
                    "edge": (r * (1 + h), r, s + 2)}[kind]
        assert degree(result) == (expected[0],) * n + (expected[1],) * n + (expected[2],) * (len(indices) * h)
    return result


KINDS = ("corona", "neighborhood", "edge")


def rook_and_shrikhande():
    vertices = list(product(range(4), repeat=2))
    rook = graph(16, [(i, j) for i, j in combinations(range(16), 2)
                      if vertices[i][0] == vertices[j][0] or vertices[i][1] == vertices[j][1]])
    steps = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    shrikhande = graph(16, [(i, j) for i, j in combinations(range(16), 2)
                           if tuple((vertices[j][k] - vertices[i][k]) % 4 for k in range(2)) in steps])
    x = sp.Symbol("x")
    expected = tuple(map(int, sp.Poly((x - 6) * (x - 2) ** 6 * (x + 2) ** 9, x).all_coeffs()))
    assert integer_charpoly(rook) == integer_charpoly(shrikhande) == expected
    assert degree(rook) == degree(shrikhande) == (6,) * 16
    k4 = lambda a: sum(all(a[i][j] for i, j in combinations(v, 2)) for v in combinations(range(16), 4))
    assert k4(rook) == 8 and k4(shrikhande) == 0  # Certifies nonisomorphism.
    print("INPUT rook/Shrikhande: degree=6, A=(x-6)(x-2)^6(x+2)^9; K4 counts=8,0")
    return rook, shrikhande


def seeds_ten():
    # Mathematical neighborhood table in Sason arXiv:2609.10884, not program code.
    rows = ["bdfh", "acde", "begj", "abfh", "bcgj", "adij", "cehi", "adgi", "fghj", "cefi"]
    a = graph(10, [(i, ord(j) - ord("a")) for i, row in enumerate(rows) for j in row])
    es = {(i, j) for i, j in combinations(range(10), 2) if a[i][j]}
    for edge in [(0, 5), (2, 6), (6, 7), (5, 9)]:
        es.remove(edge)
    es.update([(0, 6), (2, 5), (5, 7), (6, 9)])
    b = graph(10, es)
    assert degree(a) == degree(b) == (4,) * 10
    assert integer_charpoly(a) == integer_charpoly(b)
    triangles = lambda g: sorted(sum(g[i][v] * g[j][v] for v in range(10))
                                 for i, j in combinations(range(10), 2) if g[i][j])
    assert max(triangles(a)) == 3 and max(triangles(b)) == 2
    return a, b


def small_cases():
    groups = defaultdict(list)
    count = 0
    for n in range(5):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            a = graph(n, [edge for i, edge in enumerate(possible) if (mask >> i) & 1])
            d = degree(a)
            if len(set(d)) <= 1:
                groups[(n, d[0] if n else 0, integer_charpoly(a))].append(a)
                count += 1
    comparisons = 0
    quadruples = 0
    for gs in groups.values():
        for hs in groups.values():
            reference = {kind: charpoly(corona(gs[0], hs[0], kind), "NL") for kind in KINDS}
            for a in gs:
                for b in hs:
                    for kind in KINDS:
                        assert charpoly(corona(a, b, kind), "NL") == reference[kind]
                        comparisons += 1
            # Equality with a common reference checks every ordered quadruple.
            quadruples += len(gs) ** 2 * len(hs) ** 2
    print(f"SMALL exhaustive all regular labeled graphs on 0..4 vertices: graphs={count}, groups={len(groups)}, ordered factor quadruples={quadruples}, normalized reference comparisons={comparisons}")


def main():
    print("Exact arithmetic: Python integers/Fraction and SymPy", sp.__version__)
    k2 = graph(2, [(0, 1)])
    rook, shri = rook_and_shrikhande()
    # Reproduce source Theorem 5.18 / journal Theorem 39 for a genuine NICS attached pair.
    for kind in KINDS:
        left, right = corona(k2, rook, kind), corona(k2, shri, kind)
        for spectrum in ("A", "L", "Q"):
            assert charpoly(left, spectrum) == charpoly(right, spectrum)
        assert charpoly(left, "NL") == charpoly(right, "NL")
        print(f"FIDELITY Theorem 39 and expected NL: G=G'=K2, H=rook, H'=Shrikhande, {kind}, order={len(left)}, spectra=A,L,Q,NL equal")
    for kind in KINDS:
        left, right = corona(rook, k2, kind), corona(shri, k2, kind)
        assert charpoly(left, "NL") == charpoly(right, "NL")
        print(f"BASE VARIATION rook/Shrikhande with H=H'=K2: {kind}, order={len(left)}, NL charpolys equal")
    a, b = seeds_ten()
    for kind in KINDS:
        left, right = corona(a, a, kind), corona(b, b, kind)
        assert charpoly(left, "NL") == charpoly(right, "NL")
        print(f"SIMULTANEOUS VARIATION nonisomorphic 10-vertex factors: {kind}, order={len(left)}, NL charpolys equal")
    # Direct algebraic-square-root cross-check of the rational surrogate, including isolates.
    a = corona(graph(2, []), k2, "corona")
    d = degree(a)
    f = sp.diag(*[1 / sp.sqrt(v) if v else sp.Integer(0) for v in d])
    lap = sp.diag(*d) - sp.Matrix(a)
    actual = tuple(Fraction(c) for c in (f * lap * f).charpoly().all_coeffs())
    assert actual == charpoly(a, "NL") and actual[-1] == 0
    print("TOTALIZATION direct f(D)(D-A)f(D) agrees with rational surrogate including isolated duplicate vertices")
    small_cases()
    print("Finite diagnostics only; no proof of the universally quantified assertion is claimed.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
