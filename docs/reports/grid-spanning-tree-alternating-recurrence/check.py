#!/usr/bin/env python3
"""Exact controls for Raff, arXiv:0809.2551v1, Conjecture 2.

Python 3 and SymPy. No floating-point evidence. This is a finite check, not a
proof of the universal conjecture. Graphs are finite, simple, and labelled.
"""
from functools import lru_cache
from itertools import combinations
import sympy as sp

X, Y = sp.symbols("X Y")
COUNTS = {"cofactors": 0, "formula_comparisons": 0, "recurrence_equations": 0}


def laplacian(k, edges):
    a = [[0] * k for _ in range(k)]
    for u, v in edges:
        a[u][u] += 1
        a[v][v] += 1
        a[u][v] -= 1
        a[v][u] -= 1
    return a


def bareiss(a):
    """Integer determinant, with exact-division assertions and row pivoting."""
    a = [row[:] for row in a]
    size, previous, sign = len(a), 1, 1
    if not size:
        return 1
    for c in range(size - 1):
        pivot_row = next((r for r in range(c, size) if a[r][c]), None)
        if pivot_row is None:
            return 0
        if pivot_row != c:
            a[c], a[pivot_row] = a[pivot_row], a[c]
            sign = -sign
        pivot = a[c][c]
        for i in range(c + 1, size):
            for j in range(c + 1, size):
                value = pivot * a[i][j] - a[i][c] * a[c][j]
                assert value % previous == 0
                a[i][j] = value // previous
            a[i][c] = 0
        previous = pivot
    return sign * a[-1][-1]


def tree_count(k, edges, n):
    assert n >= 1
    COUNTS["cofactors"] += 1
    if k == 0:
        return 0
    product_edges = [(n * u + j, n * v + j)
                     for u, v in edges for j in range(n)]
    product_edges += [(n * u + j, n * u + j + 1)
                      for u in range(k) for j in range(n - 1)]
    a = laplacian(k * n, product_edges)
    return bareiss([row[:-1] for row in a[:-1]])


def connected(k, edges):
    if not k:
        return False
    reached = {0}
    while True:
        more = reached | {v for u, v in edges if u in reached}
        more |= {u for u, v in edges if v in reached}
        if more == reached:
            return len(reached) == k
        reached = more


@lru_cache(None)
def second_kind(n):
    """S_0=0, S_1=1, S_(n+1)=(lambda+2)S_n-S_(n-1)."""
    if n == 0:
        return sp.Integer(0)
    if n == 1:
        return sp.Integer(1)
    return sp.expand((X + 2) * second_kind(n - 1) - second_kind(n - 2))


@lru_cache(None)
def spectrum_controls(eigenvalues, tau):
    """No radical choices: resultants multiply all reciprocal-root pairs."""
    polynomial = X - 1
    for lam in eigenvalues:
        polynomial = sp.resultant(polynomial.subs(X, Y),
                                  X**2 - (lam + 2) * X * Y + Y**2, Y)
        polynomial = sp.Poly.from_list(
            [sp.simplify(c) for c in sp.Poly(polynomial, X).all_coeffs()], X).as_expr()
    polynomial = sp.Poly(polynomial, X, domain=sp.ZZ)
    coefficients = tuple(map(int, polynomial.all_coeffs()))
    assert coefficients[0] == 1
    assert polynomial.degree() == 2 ** len(eigenvalues)
    assert all((-1)**j * c > 0 for j, c in enumerate(coefficients))
    # Sturm root counting is exact and handles repeated roots by squarefree part.
    squarefree = polynomial.sqf_part()
    assert squarefree.eval(0) != 0
    assert squarefree.count_roots(0, sp.oo) == squarefree.degree()
    sequence = [0]
    for n in range(1, 20):
        value = sp.simplify(tau * sp.prod(second_kind(n).subs(X, l)
                                        for l in eigenvalues))
        assert value.is_Integer and not value.has(sp.Float)
        sequence.append(int(value))
    r = polynomial.degree()
    equations = 0
    for n in range(1, 20 - r):
        assert sum(c * sequence[n + r - j]
                   for j, c in enumerate(coefficients)) == 0
        equations += 1
    return coefficients, tuple(sequence), equations


def control(k, edges):
    if not connected(k, edges):
        return (1, -1), (0,) * 20, 18
    eig = sp.Matrix(laplacian(k, edges)).eigenvals()
    assert eig.pop(sp.Integer(0)) == 1
    eigenvalues = tuple(sorted((l for l, m in eig.items() for _ in range(m)), key=str))
    tau = bareiss([row[:-1] for row in laplacian(k, edges)[:-1]])
    return spectrum_controls(eigenvalues, tau)


def main():
    # Fidelity: Raff Theorem 6.3, with its stated T_1=1 and T_2=4.
    ladder_edges = [(0, 1)]
    ladder = [tree_count(2, ladder_edges, n) for n in range(1, 20)]
    assert ladder[:5] == [1, 4, 15, 56, 209]
    assert all(ladder[n + 2] == 4 * ladder[n + 1] - ladder[n]
               for n in range(len(ladder) - 2))
    print("FIDELITY Theorem 6.3: T_1..T_5 =", ladder[:5], "; recurrence verified through n=19")
    # General small-case family: every labelled simple graph on 0..4 vertices.
    total = connected_total = 0
    for k in range(5):
        possible_edges = list(combinations(range(k), 2))
        for mask in range(1 << len(possible_edges)):
            edges = [e for j, e in enumerate(possible_edges) if mask & (1 << j)]
            coefficients, sequence, equations = control(k, edges)
            total += 1
            connected_total += connected(k, edges)
            COUNTS["recurrence_equations"] += equations
            for n in range(1, 7):
                assert tree_count(k, edges, n) == sequence[n]
                COUNTS["formula_comparisons"] += 1
    assert (total, connected_total) == (76, 44)
    print("SMALL FAMILY: 76 labelled graphs on 0..4 vertices; 44 connected, 32 disconnected/empty")
    print("Each graph: six independent cofactors; exact spectral formula and alternating integer recurrence through n=19")
    # Source table's P_3 and P_4 initial values; G_4 recurrence coefficients.
    p3 = [(0, 1), (1, 2)]
    p4 = [(0, 1), (1, 2), (2, 3)]
    star = [(0, 1), (0, 2), (0, 3)]
    assert control(3, p3)[1][1:6] == (1, 15, 192, 2415, 30305)
    assert control(4, p4)[1][1:6] == (1, 56, 2415, 100352, 4140081)
    expected = {
        "P4": (1, -56, 672, -2632, 4094, -2632, 672, -56, 1),
        "K1,3": (1, -54, 625, -2646, 4164, -2646, 625, -54, 1),
    }
    for name, edges in [("P4", p4), ("K1,3", star)]:
        coefficients, sequence, _ = control(4, edges)
        assert coefficients == expected[name]
        for n in range(7, 20):
            assert tree_count(4, edges, n) == sequence[n]
            COUNTS["formula_comparisons"] += 1
        print(name, "annihilator descending coefficients:", list(coefficients))
    # Independent complete-graph controls, including larger repeated spectra.
    for k in range(2, 9):
        edges = list(combinations(range(k), 2))
        for n in range(1, 7):
            sn = second_kind(n).subs(X, k)
            value = k ** (k - 2) * sn ** (k - 1)
            assert not value.has(sp.Float)
            assert tree_count(k, edges, n) == value
            COUNTS["formula_comparisons"] += 1
    print("COMPLETE FAMILY: 42 cofactors for K_k box P_n, 2<=k<=8, 1<=n<=6")
    print("Exact arithmetic: Python integers and SymPy", sp.__version__)
    print("COUNTS", COUNTS)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
