#!/usr/bin/env python3
"""Exact checks for arXiv:2112.09744v1, Conjecture 4.1.

Requires SymPy 1.14.0. No numerical root approximations are used.
"""
from itertools import combinations
from math import comb

import sympy as sp


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def graph_from_edges(n, edges):
    adjacency = [0] * n
    for a, b in edges:
        require(0 <= a < b < n, "simple undirected edge")
        adjacency[a] |= 1 << b
        adjacency[b] |= 1 << a
    return adjacency


def connected(adjacency, deleted=0):
    remaining = ((1 << len(adjacency)) - 1) & ~deleted
    if not remaining:
        return False
    frontier = remaining & -remaining
    reached = 0
    while frontier:
        reached |= frontier
        next_frontier = 0
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            next_frontier |= adjacency[bit.bit_length() - 1]
        frontier = next_frontier & remaining & ~reached
    return reached == remaining


def clique_counts(adjacency):
    """Enumerate every clique once by extending its increasing vertex list."""
    counts = [1]

    def visit(candidates, size):
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            if len(counts) <= size:
                counts.append(0)
            counts[size] += 1
            visit(candidates & adjacency[bit.bit_length() - 1], size + 1)

    visit((1 << len(adjacency)) - 1, 1)
    return counts


def triangle_free(adjacency):
    for a, neighbors in enumerate(adjacency):
        for b in range(a + 1, len(adjacency)):
            if (neighbors >> b) & 1 and neighbors & adjacency[b]:
                return False
    return True


def fidelity():
    """Reproduce Proposition 3.3 on all labelled graphs with at most 6 vertices."""
    print("FIDELITY: source Proposition 3.3 (connected triangle-free graphs)")
    for n in range(1, 7):
        pairs = list(combinations(range(n), 2))
        checked = 0
        for mask in range(1 << len(pairs)):
            adjacency = graph_from_edges(n, (e for i, e in enumerate(pairs) if mask >> i & 1))
            if not triangle_free(adjacency) or not connected(adjacency):
                continue
            counts = clique_counts(adjacency)
            edges = sum(neighbors.bit_count() for neighbors in adjacency) // 2
            require(counts == ([1, n] if edges == 0 else [1, n, edges]), "triangle-free counts")
            if edges:
                require(n * n - 4 * edges >= 0, "quadratic has only real roots")
            checked += 1
        require(checked > 0, "nonempty fidelity set")
        print(f"n={n}: all {checked} connected triangle-free labelled graphs pass")
    x = sp.Symbol("x")
    for n in range(2, 11):
        adjacency = graph_from_edges(n, ((i, i + 1) for i in range(n - 1)))
        p = sum(c * x**k for k, c in enumerate(clique_counts(adjacency)))
        require(sp.expand(p - (1 + x) * (1 + (n - 1) * x)) == 0, "source tree factorization")
        require(p.subs(x, -1) == p.subs(x, -sp.Rational(1, n - 1)) == 0, "source tree roots")
    print("Source Proposition 3.1: paths on 2..10 vertices have roots -1 and -1/(n-1)")


def counterfamily(m):
    n = 6 * m + 1
    core_size = 6 * m
    return graph_from_edges(n, (
        (a, b) for a, b in combinations(range(n), 2)
        if (b < core_size and a // (2 * m) != b // (2 * m))
        or (b == core_size and a < m)
    ))


def family_checks():
    x, u = sp.symbols("x u")
    m = sp.Symbol("m", integer=True, positive=True)
    p = (1 + 2 * m * x)**3 + x + m * x**2
    expected = 1 + (6 * m + 1) * x + (12 * m**2 + m) * x**2 + 8 * m**3 * x**3
    require(sp.expand(p - expected) == 0, "uniform coefficient identity")
    discriminant = sp.factor(sp.discriminant(p, x))
    require(sp.expand(discriminant - m**2 * (1 - 108 * m**2)) == 0, "uniform discriminant identity")
    q = 4 * m * u**3 + u**2 - 1
    require(sp.cancel(4 * m * p.subs(x, (u - 1) / (2 * m)) - q) == 0, "uniform affine identity")
    require(sp.expand(sp.diff(q, u) - 2 * u * (6 * m * u + 1)) == 0, "derivative")
    require(sp.cancel(q.subs(u, -1 / (6 * m)) - (1 / (108 * m**2) - 1)) == 0, "critical value")
    require(q.subs(u, 0) == -1, "second critical value")
    # Exact coefficient sign proves 108*m^2-1 >= 431 for all integers m>=2.
    t = sp.Symbol("t", integer=True, nonnegative=True)
    require(sp.Poly((108 * m**2 - 1).subs(m, t + 2) - 431, t).all_coeffs() == [108, 432, 0], "uniform negative discriminant for m>=2")
    print("FAMILY: symbolic identities in integer parameter m")
    print("C_m(x) = 1 + (6*m+1)*x + (12*m**2+m)*x**2 + 8*m**3*x**3")
    print(f"discriminant = {discriminant}")
    print("4*m*C_m((u-1)/(2*m)) = 4*m*u**3 + u**2 - 1")
    print("critical values: 1/(108*m**2)-1 and -1; both negative for m>=2")
    for parameter in range(2, 9):
        adjacency = counterfamily(parameter)
        counts = clique_counts(adjacency)
        require(counts == [1, 6 * parameter + 1, 12 * parameter**2 + parameter, 8 * parameter**3], "adjacency-based complete clique enumeration")
        # Independent direct 4-subset check, separate from recursive enumeration.
        require(not any(all(adjacency[a] >> b & 1 for a, b in combinations(s, 2))
                        for s in combinations(range(len(adjacency)), 4)), "K4-free")
        polynomial = sp.Poly(p.subs(m, parameter), x, domain=sp.QQ)
        require(polynomial.count_roots(-sp.oo, sp.oo) == 1, "exact Sturm real-root count")
        require(sp.gcd(polynomial, polynomial.diff()).degree() == 0, "all three roots are simple")
        require(polynomial.degree() == 3, "one real plus two nonreal conjugate roots")
        cut = (1 << parameter) - 1
        require(cut.bit_count() == parameter and not connected(adjacency, cut), "upper bound on connectivity")
        deletion_checks = 0
        if parameter <= 5:
            for k in range(parameter):
                for subset in combinations(range(len(adjacency)), k):
                    deleted = sum(1 << v for v in subset)
                    require(connected(adjacency, deleted), "survives every deletion of fewer than m vertices")
                    deletion_checks += 1
            require(deletion_checks == sum(comb(len(adjacency), k) for k in range(parameter)), "exhaustive deletion count")
        print(f"m={parameter}: clique counts={counts}; K4-free; exactly 1 real root; cut size={parameter}; deletion checks={deletion_checks}")
    print("Finite graph checks cover m=2..8; exhaustive lower connectivity checks cover m=2..5.")
    print("Uniform graph counting/connectivity and root conclusions still require mathematical proofs.")


def main():
    print(f"Exact arithmetic; SymPy {sp.__version__}; no floating-point evidence")
    fidelity()
    family_checks()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
