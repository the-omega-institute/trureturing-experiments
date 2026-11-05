#!/usr/bin/env python3
"""Exact recurrence, source-fidelity, and finite family checks; Python stdlib only."""
from collections import deque
from fractions import Fraction
from itertools import combinations


def creation(alpha):
    """Source recurrence: complement of (old graph disjoint-union a-clique)."""
    assert alpha and all(a >= 1 for a in alpha)
    graph, blocks = [], []
    for a in alpha:
        n = len(graph)
        union = [[0] * (n + a) for _ in range(n + a)]
        for i in range(n):
            union[i][:n] = graph[i][:]
        for i in range(n, n + a):
            for j in range(n, n + a):
                union[i][j] = int(i != j)
        graph = [[int(i != j and not union[i][j])
                  for j in range(n + a)] for i in range(n + a)]
        blocks.append(list(range(n, n + a)))
    return graph, blocks


def distance_signless(graph):
    """BFS distances, independent of the proposed diameter-two shortcut."""
    n = len(graph)
    distances = []
    for start in range(n):
        d, todo = [None] * n, deque([start])
        d[start] = 0
        while todo:
            u = todo.popleft()
            for v in range(n):
                if graph[u][v] and d[v] is None:
                    d[v] = d[u] + 1
                    todo.append(v)
        assert None not in d, 'disconnected source graph'
        distances.append(d)
    transmissions = [sum(row) for row in distances]
    q = [[distances[i][j] + (transmissions[i] if i == j else 0)
          for j in range(n)] for i in range(n)]
    return q, transmissions, distances


def mv(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    pivot = 0
    for col in range(len(a[0])):
        row = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if row is None:
            continue
        a[pivot], a[row] = a[row], a[pivot]
        lead = a[pivot][col]
        a[pivot] = [x / lead for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            lead = a[i][col]
            a[i] = [x - lead * y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def charpoly(a):
    """Faddeev--LeVerrier over Z, checking every exact division; descending coefficients."""
    n = len(a)
    b = [[int(i == j) for j in range(n)] for i in range(n)]
    coefficients = [1]
    for k in range(1, n + 1):
        ab = [[sum(a[i][h] * b[h][j] for h in range(n))
               for j in range(n)] for i in range(n)]
        numerator = -sum(ab[i][i] for i in range(n))
        assert numerator % k == 0
        c = numerator // k
        coefficients.append(c)
        b = [[ab[i][j] + (c if i == j else 0) for j in range(n)]
             for i in range(n)]
    return coefficients


def poly_product(*polys):
    result = [1]
    for p in polys:
        out = [0] * (len(result) + len(p) - 1)
        for i, a in enumerate(result):
            for j, b in enumerate(p):
                out[i + j] += a * b
        result = out
    return result


def vicinal(graph, u, v):
    return all(not edge or w == v or graph[v][w]
               for w, edge in enumerate(graph[u]))


def is_chain(graph, vertices):
    return all(vicinal(graph, u, v) or vicinal(graph, v, u)
               for u, v in combinations(vertices, 2))


def exceptions(alpha, blocks, transmission):
    result = set()
    for i, (a, block) in enumerate(zip(alpha, blocks), 1):
        assert len({transmission[v] for v in block}) == 1
        if a >= 2:
            result.add(transmission[block[0]] - (1 if i % 2 else 2))
    return result


def fidelity():
    alpha = [1, 3, 1, 2]
    graph, blocks = creation(alpha)
    q, t, _ = distance_signless(graph)
    quotient = [[sum(q[block[0]][v] for v in target)
                 for target in blocks] for block in blocks]
    expected = [[7, 3, 2, 2], [1, 13, 2, 2],
                [2, 6, 10, 2], [1, 3, 1, 9]]
    assert quotient == expected
    for block, row in zip(blocks, quotient):
        for u in block:
            assert [sum(q[u][v] for v in b) for b in blocks] == row
    assert charpoly(quotient) == poly_product([1, -8], [1, -8], [1, -23, 100])
    assert charpoly(q) == poly_product([1, -7], [1, -7], [1, -5],
                                      [1, -8], [1, -8], [1, -23, 100])
    assert 8 not in exceptions(alpha, blocks, t)
    all_vertices = set(range(len(graph)))
    covers = [(sorted(s), sorted(all_vertices - s))
              for size in range(1, len(graph))
              for subset in combinations(range(len(graph)), size)
              for s in [set(subset)]
              if is_chain(graph, s) and is_chain(graph, all_vertices - s)]
    assert covers and not is_chain(graph, all_vertices)
    incomparable = next((u, v) for u, v in combinations(range(len(graph)), 2)
                        if not vicinal(graph, u, v) and not vicinal(graph, v, u))
    print('FIDELITY Example 5.4: quotient =', quotient)
    print('FIDELITY quotient charpoly = (L-8)^2*(L^2-23*L+100)')
    print('FIDELITY eigenvalue 8 has algebraic multiplicity 2; TQ =', sorted(exceptions(alpha, blocks, t)))
    print('FIDELITY Dilworth = 2; 0-based chain cover =', covers[0],
          '; incomparable pair =', incomparable)


def family(r):
    R, k, n, lam = r + 1, r + 2, 3 * r + 5, 4 * r + 4
    alpha = [r] + [1] * (2 * r + 5)
    graph, blocks = creation(alpha)
    x = list(range(R))
    p = [R + 2 * j for j in range(k)]
    q = [R + 2 * j + 1 for j in range(k)]
    # Independently reconstruct the isolated/dominating description.
    direct = [[0] * n for _ in range(n)]
    clique = x + q
    for a, b in combinations(clique, 2):
        direct[a][b] = direct[b][a] = 1
    for j in range(k):
        for h in range(j, k):
            direct[p[j]][q[h]] = direct[q[h]][p[j]] = 1
    assert graph == direct
    matrix, t, distances = distance_signless(graph)
    assert all(graph[q[-1]][v] for v in range(n) if v != q[-1])
    assert all(distances[i][j] == (0 if i == j else 1 if graph[i][j] else 2)
               for i in range(n) for j in range(n))
    assert [t[v] for v in x] == [4 * r + 6] * R
    assert [t[v] for v in p] == [5 * r + 5 + j for j in range(1, k + 1)]
    assert [t[v] for v in q] == [4 * r + 6 - j for j in range(1, k + 1)]
    assert t == [2 * (n - 1) - sum(row) for row in graph]
    u, v = [0] * n, [0] * n
    for a in x:
        u[a] = -1
    u[p[0]], u[q[0]], v[q[0]] = 1, r, R
    for j in range(1, k):
        v[p[j]], v[q[j]] = -1, 1
    assert mv(matrix, u) == [lam * z for z in u]
    assert mv(matrix, v) == [lam * z for z in v]
    minor = u[x[0]] * v[p[1]] - u[p[1]] * v[x[0]]
    assert minor == 1
    shifted = [[matrix[i][j] - (lam if i == j else 0)
                for j in range(n)] for i in range(n)]
    nullity = n - rank(shifted)
    assert nullity >= 2
    chain = list(reversed(p)) + x + q
    assert sorted(chain) == list(range(n))
    assert all(vicinal(graph, chain[i], chain[j])
               for i in range(n) for j in range(i, n))
    tq = exceptions(alpha, blocks, t)
    assert tq == (set() if r == 1 else {4 * r + 5}) and lam not in tq
    # The maximal true-twin class also has a different local eigenvalue.
    assert t[x[0]] - 1 == lam + 1
    print(f'FAMILY r={r}: order={n}, lambda={lam}, nullity={nullity}, '
          f'Dilworth=1, TQ={sorted(tq)}, independence_minor={minor}; PASS')
    if r == 1:
        expected = poly_product([1, -9], [1, -8], [1, -8],
                                [1, -55, 1139, -11261, 53580, -98676])
        assert charpoly(matrix) == expected
        assert sum(map(sum, graph)) // 2 == 16
        assert t == [10, 10, 11, 9, 12, 8, 13, 7]
        print('BASE r=1: Q =', matrix)
        print('BASE r=1: u =', u, '; v =', v)
        print('BASE r=1: Qu =', mv(matrix, u), '; Qv =', mv(matrix, v))
        print('BASE r=1: 1-based vicinal chain =', [i + 1 for i in chain])
        print('BASE r=1: charpoly = (L-9)*(L-8)^2*'
              '(L^5-55*L^4+1139*L^3-11261*L^2+53580*L-98676)')


if __name__ == '__main__':
    fidelity()
    for parameter in range(1, 13):
        family(parameter)
    print('Scope: exact finite checks r=1,...,12; the infinite family still requires proof.')
    print('ALL_CHECKS_OK')
