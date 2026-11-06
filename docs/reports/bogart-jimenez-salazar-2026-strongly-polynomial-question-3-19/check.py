#!/usr/bin/env python3
"""Exact independent checks for Bogart--Jimenez-Salazar Question 3.19.

Only Python's standard library is used. All evidence is integer arithmetic.
Finite checks support the preregistered argument; they do not prove its
unbounded quantifiers or certify literature novelty.
"""
from itertools import combinations, product
from math import comb


def multipartite(sizes):
    parts = tuple(i for i, size in enumerate(sizes) for _ in range(size))
    return tuple(tuple(a != b for b in parts) for a in parts)


def brute_hom(source, target):
    edges = [(i, j) for i in range(len(source)) for j in range(i)
             if source[i][j]]
    return sum(all(target[f[i]][f[j]] for i, j in edges)
               for f in product(range(len(target)), repeat=len(source)))


def common_neighbor_count(m, target):
    return sum(sum(target[u][x] and target[v][x]
                   for x in range(len(target))) ** m
               for u in range(len(target)) for v in range(len(target)))


def proposed_count(m, n):
    return 2 * n**2 * (n + 1)**m + 4 * n**(m + 1) + 2**m * n**m + 2 * n**2


def coefficients(m):
    """Expand the proposed expression without interpolation or floating point."""
    c = [0] * (m + 3)
    for j in range(m + 1):
        c[j + 2] += 2 * comb(m, j)
    c[m + 1] += 4
    c[m] += 2**m
    c[2] += 2
    return c


def coloring_sum(source, sizes):
    edges = [(i, j) for i in range(len(source)) for j in range(i)
             if source[i][j]]
    total = 0
    for coloring in product(range(len(sizes)), repeat=len(source)):
        if all(coloring[i] != coloring[j] for i, j in edges):
            lifts = 1
            for color in coloring:
                lifts *= sizes[color]
            total += lifts
    return total


def all_graphs(v):
    pairs = list(combinations(range(v), 2))
    for mask in range(1 << len(pairs)):
        a = [[False] * v for _ in range(v)]
        for bit, (i, j) in enumerate(pairs):
            a[i][j] = a[j][i] = bool(mask & (1 << bit))
        yield tuple(map(tuple, a))


def main():
    # Fidelity: Proposition 3.7(1), stars on q+1 and n+1 vertices.
    fidelity = 0
    for q in range(6):
        for n in range(5):
            assert brute_hom(multipartite((1, q)), multipartite((1, n))) == n**q + n
            fidelity += 1
    print(f'FIDELITY Proposition 3.7(1): {fidelity} exact full-map enumerations OK')

    enumerations = 0
    for m in range(3, 6):
        for n in range(3):
            actual = brute_hom(multipartite((m, 2)), multipartite((1, n, n)))
            assert actual == proposed_count(m, n), (m, n, actual)
            enumerations += 1
            print(f'FULL_MAP m={m} n={n} hom={actual}')
    print(f'CONSTRUCTION: {enumerations} full-map enumerations OK')

    checked = 0
    for m in range(3, 9):
        for n in range(7):
            assert common_neighbor_count(m, multipartite((1, n, n))) == proposed_count(m, n)
            checked += 1
    print(f'COMMON_NEIGHBORS: {checked} exact cases, 3<=m<=8 and 0<=n<=6, OK')

    graphs = evaluations = 0
    for v in range(5):
        for source in all_graphs(v):
            graphs += 1
            for parameter in range(4):
                for sizes in ((parameter + 3, 2), (1, parameter, parameter)):
                    assert brute_hom(source, multipartite(sizes)) == coloring_sum(source, sizes)
                    evaluations += 1
    print(f'COLORING_FIBERS: {graphs} labelled graphs on <=4 vertices; {evaluations} evaluations OK')

    for q in range(31):
        m = q + 3
        c = coefficients(m)
        assert len(c) - 1 == q + 5 and c[-1] == 2
        assert c[-2] == 2*q + 10
        assert c[-3] == (q + 3)*(q + 2) + 2**(q + 3)
        for n in range(7):
            assert sum(a * n**j for j, a in enumerate(c)) == proposed_count(m, n)
    print('COEFFICIENTS: q=0..30, degree q+5, c0=2, c1=2q+10, c2=(q+3)(q+2)+2^(q+3) OK')

    residual = [coefficients(q+3)[-3] - (q+3)*(q+2) for q in range(32)]
    for k in range(11):
        assert residual == [2**(q+3) for q in range(len(residual))]
        residual = [b-a for a, b in zip(residual, residual[1:])]
    print('FORWARD_DIFFERENCES: orders 0..10 of residual equal 2^(q+3), exactly OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
