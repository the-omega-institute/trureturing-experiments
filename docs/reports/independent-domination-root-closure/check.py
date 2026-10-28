#!/usr/bin/env python3
"""Exact finite checks for Shao--Zhao, arXiv:2609.06569v5, Open Problem 1.

Standard-library integers and fractions only. No numerical root approximation.
The finite checks do not prove the two topological closure statements.
"""
from fractions import Fraction
from itertools import combinations


def graph(n, edges):
    """An undirected simple graph as a tuple of integer neighborhood masks."""
    neighbors = [0] * n
    for u, v in edges:
        assert 0 <= u < v < n
        neighbors[u] |= 1 << v
        neighbors[v] |= 1 << u
    return tuple(neighbors)


def polynomial(neighbors):
    """Enumerate actual subsets; test independence and outside-neighbor domination."""
    n = len(neighbors)
    full = (1 << n) - 1
    counts = [0] * (n + 1)
    for subset in range(1 << n):
        covered = subset
        independent = True
        for v, adjacent in enumerate(neighbors):
            if subset & (1 << v):
                if adjacent & subset:
                    independent = False
                    break
                covered |= adjacent
        if independent and covered == full:
            counts[subset.bit_count()] += 1
    return counts


def clique_substitution(neighbors, r):
    n = len(neighbors)
    return graph(n * r, ((u, v) for u, v in combinations(range(n * r), 2)
                        if u // r == v // r or
                        (neighbors[u // r] & (1 << (v // r)))))


def multipartite(parts):
    labels = tuple(i for i, size in enumerate(parts) for _ in range(size))
    return graph(len(labels), ((u, v) for u, v in combinations(range(len(labels)), 2)
                               if labels[u] != labels[v]))


def partitions(n, lower=1):
    if n == 0:
        yield ()
    for first in range(lower, n + 1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def evaluate(coefficients, value):
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def connected(neighbors):
    if not neighbors:
        return False
    reached = 1
    while True:
        grown = reached
        for v, adjacent in enumerate(neighbors):
            if reached & (1 << v):
                grown |= adjacent
        if grown == reached:
            return reached == (1 << len(neighbors)) - 1
        reached = grown


def check_fidelity():
    """Source Section 2 (also Open Problem 1): D_i(G[K_r],x) = D_i(G,r*x)."""
    cases = 0
    for n in range(5):
        edges = tuple(combinations(range(n), 2))
        for mask in range(1 << len(edges)):
            g = graph(n, (edge for i, edge in enumerate(edges) if mask & (1 << i)))
            base = polynomial(g)
            for r in range(1, 4):
                actual = polynomial(clique_substitution(g, r))
                expected = [0] * (n * r + 1)
                for k, coefficient in enumerate(base):
                    expected[k] = coefficient * r ** k
                assert actual == expected, (n, mask, r, actual, expected)
                cases += 1
    print(f'FIDELITY scaling identity: {cases} cases; every labelled G on 0..4 vertices, r=1..3')


def check_multipartite():
    cases = 0
    for n in range(1, 13):
        for parts in partitions(n):
            actual = polynomial(multipartite(parts))
            expected = [0] * (n + 1)
            for size in parts:
                expected[size] += 1
            assert actual == expected, (parts, actual, expected)
            cases += 1
    print(f'MULTIPARTITE identity: {cases} cases; all positive part-size partitions of 1..12')


def check_witness_family():
    cases = 0
    real_cases = 0
    for m in range(1, 6):
        for a in range(1, 15):
            for b in range(1, 8):
                n = m * (a + 2 * b)
                if n > 14:
                    continue
                g = multipartite((m,) * a + (2 * m,) * b)
                assert n >= 3 and connected(g)
                actual = polynomial(g)
                expected = [0] * (n + 1)
                expected[m] = a
                expected[2 * m] = b
                assert actual == expected, (a, b, m, actual, expected)
                assert evaluate(actual, Fraction(0)) == 0
                # Exact reduction modulo b*x^m+a, without extracting complex roots.
                remainder = [Fraction(c) for c in actual]
                for degree in range(n, m - 1, -1):
                    quotient = remainder[degree] / b
                    remainder[degree] = Fraction(0)
                    remainder[degree - m] -= quotient * a
                assert all(c == 0 for c in remainder)
                if m == 1:
                    root = -Fraction(a, b)
                    assert root < 0 and evaluate(actual, root) == 0
                    real_cases += 1
                cases += 1
    print(f'WITNESS identity a*x^m+b*x^(2m): {cases} connected graphs, m*(a+2*b)<=14')
    print(f'REAL roots -a/b: {real_cases} exact Fraction evaluations')
    # An exact nonreal example: D_i(K_{2,4},x)=x^2+x^4 and i^2=-1.
    p = polynomial(multipartite((2, 4)))
    real = sum(c * (-1) ** (k // 2) for k, c in enumerate(p) if k % 2 == 0)
    imag = sum(c * (-1) ** ((k - 1) // 2) for k, c in enumerate(p) if k % 2 == 1)
    assert (real, imag) == (0, 0)
    print('COMPLEX example K_{2,4}: evaluation at i is exactly (0,0) in Z[i]')


def check_positive_and_empty():
    assert polynomial(()) == [1]
    cases = 0
    for n in range(1, 6):
        edges = tuple(combinations(range(n), 2))
        for mask in range(1 << len(edges)):
            p = polynomial(graph(n, (e for i, e in enumerate(edges) if mask & (1 << i))))
            assert p[0] == 0 and sum(p) > 0 and all(c >= 0 for c in p)
            for value in (Fraction(1, 7), Fraction(1), Fraction(5, 2)):
                assert evaluate(p, value) > 0
            cases += 1
    print(f'POSITIVE checks: {cases} nonempty labelled graphs on 1..5 vertices, x=1/7,1,5/2')
    print('EMPTY graph: polynomial 1, hence no roots')


def main():
    check_fidelity()
    check_multipartite()
    check_witness_family()
    check_positive_and_empty()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
