#!/usr/bin/env python3
"""Exact, standard-library checks for KMV arXiv:2610.22248v1, Question 27.

Structure matrices use ROWS: table[i] is e_i squared in the old basis.
Change-of-basis matrices use COLUMNS: P[:,i] is f_i in the old basis.
No floating point, random sampling, or Lean proof is used.
"""

from collections import deque
from fractions import Fraction
from math import gcd


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def inverse_and_det(a):
    n = len(a)
    b = [[Fraction(x) for x in row] + list(map(Fraction, unit))
         for row, unit in zip(a, identity(n))]
    det = Fraction(1)
    for j in range(n):
        pivot = next(i for i in range(j, n) if b[i][j])
        if pivot != j:
            b[j], b[pivot] = b[pivot], b[j]
            det = -det
        value = b[j][j]
        det *= value
        b[j] = [x / value for x in b[j]]
        for i in range(n):
            if i != j and b[i][j]:
                value = b[i][j]
                b[i] = [x - value * y for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b], det


def product(table, x, y):
    return [sum(x[i] * y[i] * table[i][j] for i in range(len(table)))
            for j in range(len(table))]


def transport(table, p):
    """Compute every f_i f_j by bilinearity, then return to F coordinates."""
    q, det = inverse_and_det(p)
    n = len(table)
    assert matmul(q, p) == identity(n) == matmul(p, q)
    columns = list(zip(*p))
    new = []
    for i in range(n):
        for j in range(n):
            old = product(table, columns[i], columns[j])
            coeff = [sum(q[a][b] * old[b] for b in range(n)) for a in range(n)]
            if i != j:
                assert coeff == [0] * n, (i, j, coeff)
            else:
                assert all(x.denominator == 1 for x in coeff)
                new.append([int(x) for x in coeff])
    return new, det


def rank(a, p):
    """Row rank over Q (p=0) or the prime field F_p, exactly."""
    a = [[Fraction(x) for x in row] if p == 0 else [x % p for x in row] for row in a]
    row = 0
    pivots = []
    for j in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        inv = 1 / a[row][j] if p == 0 else pow(a[row][j], -1, p)
        a[row] = [x * inv if p == 0 else (x * inv) % p for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][j]:
                value = a[i][j]
                a[i] = [x - value * y if p == 0 else (x - value * y) % p
                        for x, y in zip(a[i], a[row])]
        pivots.append(j)
        row += 1
    return row, pivots


def annihilator_check(table, p):
    n = len(table)
    # Coefficients of x*e_j, for every e_j and every output coordinate.
    equations = [[table[j][k] if i == j else 0 for i in range(n)]
                 for j in range(n) for k in range(n)]
    r, pivots = rank(equations, p)
    assert r == n - 1 and pivots == list(range(n - 1))
    assert table[-1] == [0] * n


def edges(table, p):
    return {(i, j) for i, row in enumerate(table) for j, x in enumerate(row)
            if (x != 0 if p == 0 else x % p != 0)}


def adjacency(n, arcs):
    a = [[] for _ in range(n)]
    for u, v in sorted(arcs):
        # Literal displayed source formula: a loop occurs in both counts.
        weight = int(u != v)
        a[u].append((v, weight))
        if u != v:
            a[v].append((u, -weight))
    return a


def balance_by_fundamental_cycles(n, arcs):
    """gcd of ALL closed-walk balances via an integral spanning forest.

    Tree paths give potentials h. Each arc has defect h[u]+w-h[v], the
    balance of its fundamental closed walk. On every closed walk the h
    terms telescope, so its balance is an integer combination of defects.
    Thus these defects generate exactly the closed-walk balance subgroup.
    Antiparallel arrows remain distinct edges; loops have literal weight 0.
    """
    a = adjacency(n, arcs)
    h = {}
    components = 0
    for root in range(n):
        if root in h:
            continue
        components += 1
        h[root] = 0
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v, w in a[u]:
                if v not in h:
                    h[v] = h[u] + w
                    queue.append(v)
    b = 0
    for u, v in arcs:
        b = gcd(b, abs(h[u] + int(u != v) - h[v]))
    return b, components


def balance_by_bounded_walks(n, arcs, max_length):
    """Separate exact enumeration; a cross-check, not an all-length claim."""
    a = adjacency(n, arcs)
    b = 0
    for root in range(n):
        states = {(root, 0)}
        for _ in range(max_length):
            states = {(v, score + w) for u, score in states for v, w in a[u]}
            for v, score in states:
                if v == root:
                    b = gcd(b, abs(score))
    return b


def family(m):
    n = m + 2
    table = [[0] * n for _ in range(n)]
    for i in range(m - 1):
        table[i][i + 1] = 1
    table[m - 1][0] = table[m - 1][m] = table[m][m + 1] = 1
    p = identity(n)
    p[m + 1][0] = 1
    new, det = transport(table, p)
    expected = [row[:] for row in table]
    expected[m - 1][m + 1] = -1
    assert det == 1 and new == expected
    return table, new


def fidelity():
    # Survey Section 1.1, third example: e1^2=e2, e2^2=0;
    # nilpotent of index 3 and graph a single arrow.
    table = [[0, 1], [0, 0]]
    basis = identity(2)
    assert edges(table, 0) == {(0, 1)}
    a2 = {tuple(product(table, x, y)) for x in basis for y in basis}
    assert (0, 1) in a2
    assert all(product(table, list(x), y) == [0, 0] for x in a2 for y in basis)
    assert balance_by_fundamental_cycles(2, edges(table, 0)) == (0, 1)
    print('FIDELITY survey Section 1.1: single arrow, A^2 != 0, A^3 = 0, index 3')
    # Independent sanity checks for the graph routine, including the source's
    # literal loop formula and opposing arrows with distinct identities.
    assert balance_by_fundamental_cycles(1, {(0, 0)}) == (0, 1)
    assert balance_by_fundamental_cycles(2, {(0, 1), (1, 0)}) == (2, 1)
    assert balance_by_fundamental_cycles(3, {(0, 1), (1, 2), (0, 2)}) == (1, 1)
    print('GRAPH SANITY literal-loop=0, directed-2-cycle=2, shortcut-triangle=1')


def main():
    fidelity()
    checked = 0
    for m in range(2, 21):
        old, new = family(m)
        n = m + 2
        for p in (0, 2, 3, 5, 7):
            annihilator_check(old, p)
            annihilator_check(new, p)
            e_old, e_new = edges(old, p), edges(new, p)
            assert e_new == e_old | {(m - 1, m + 1)}
            assert all(u != v for u, v in e_new)
            assert balance_by_fundamental_cycles(n, e_old) == (m, 1)
            assert balance_by_fundamental_cycles(n, e_new) == (1, 1)
            heights = list(range(m)) + [0, 1]
            assert all((heights[v] - heights[u]) % m == 1 for u, v in e_old)
            assert {(m - 1, m), (m, m + 1), (m - 1, m + 1)} <= e_new
            if m <= 6:
                assert balance_by_bounded_walks(n, e_old, 7) == m
                assert balance_by_bounded_walks(n, e_new, 7) == 1
            checked += 1
        print(f'FAMILY m={m:2d}, dim={n:2d}, Ann=span(z), det=1, '
              f'connected=both, balances={m},1; fields=Q,F2,F3,F5,F7')
    old, new = family(2)
    print('EXPLICIT m=2 row matrix B =', old)
    print('EXPLICIT m=2 row matrix F =', new)
    print(f'EXACT FAMILY CHECKS {checked}; m=2..20; bounded-walk cross-check m=2..6, lengths<=7')
    print('SCOPE finite computations support the preregistration; the all-field/all-m argument is prose, not Lean-verified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
