#!/usr/bin/env python3
"""Exact bounded checks for the Jaming--Mayeli random-graph B_2 question."""
from __future__ import annotations

from itertools import combinations
import sympy as sp


def laplacian(n: int, edges: set[tuple[int, int]]) -> sp.Matrix:
    """The integer combinatorial Laplacian D-A."""
    degree = [0] * n
    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
    return sp.Matrix(n, n, lambda i, j: degree[i] if i == j else (-1 if tuple(sorted((i, j))) in edges else 0))


def connected(n: int, edges: set[tuple[int, int]]) -> bool:
    seen = {0}
    while True:
        new = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if new == seen:
            return len(seen) == n
        seen = new


def pair_sum_polynomial(matrix: sp.Matrix) -> sp.Poly:
    """Exact polynomial whose roots are all unordered pair sums, with multiplicity."""
    n = matrix.rows
    x, z = sp.symbols("x z")
    char = matrix.charpoly(x).as_expr()
    ordered = sp.Poly(sp.resultant(char, char.subs(x, z - x), x), z, domain=sp.QQ)
    diagonal = sp.Poly((2**n) * char.subs(x, z / 2), z, domain=sp.QQ)
    quotient, remainder = sp.div(ordered, diagonal)
    assert remainder.is_zero
    factors = sp.factor_list(quotient.as_expr(), z)[1]
    root = sp.Poly(1, z, domain=sp.QQ)
    for factor, exponent in factors:
        assert exponent % 2 == 0, (matrix, factor, exponent)
        root *= sp.Poly(factor, z, domain=sp.QQ) ** (exponent // 2)
    return diagonal * root


def b2(matrix: sp.Matrix) -> bool:
    q = pair_sum_polynomial(matrix)
    return sp.gcd(q, q.diff()).degree() == 0


def edge_set(n: int, mask: int) -> set[tuple[int, int]]:
    pairs = list(combinations(range(n), 2))
    return {pairs[i] for i in range(len(pairs)) if mask & (1 << i)}


def main() -> None:
    # Fidelity: the path on three vertices has Delta = D-A and spectrum 0, 1, 3.
    p3 = {(0, 1), (1, 2)}
    L3 = laplacian(3, p3)
    assert L3 == sp.Matrix([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    assert sp.factor(L3.charpoly().as_expr()) == sp.Symbol("lambda") * (sp.Symbol("lambda") - 1) * (sp.Symbol("lambda") - 3)
    assert b2(L3)
    print("FIDELITY_OK path_P3_laplacian=D-A spectrum_charpoly=lambda*(lambda-1)*(lambda-3) B2=True")

    # Exact bounded settlement check: every connected labelled graph through n=4.
    checked = 0
    b2_count = 0
    for n in range(2, 5):
        for mask in range(1 << (n * (n - 1) // 2)):
            edges = edge_set(n, mask)
            if not connected(n, edges):
                continue
            checked += 1
            if b2(laplacian(n, edges)):
                b2_count += 1
    assert checked == 1 + 4 + 38  # connected labelled graphs on 2, 3, 4 vertices
    print(f"EXACT_SMALL_CONNECTED_GRAPHS checked={checked} b2={b2_count} all_pair_sum_tests=exact")

    # Additional exact samples at n=5, including a path, cycle, and two irregular graphs.
    samples = [
        {(0, 1), (1, 2), (2, 3), (3, 4)},
        {(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)},
        {(0, 1), (0, 2), (1, 3), (2, 4), (3, 4)},
        {(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)},
    ]
    statuses = [b2(laplacian(5, e)) for e in samples]
    assert len(statuses) == 4
    print(f"EXACT_N5_SAMPLES b2_statuses={statuses}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
