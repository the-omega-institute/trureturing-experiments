#!/usr/bin/env python3
"""Exact small-case tests for Olteanu, arXiv:2003.10187v2, Question 5.2.

Compute projective dimension of S/I(G) using Hochster's formula and reduced
simplicial boundary ranks. Use Q, F_2 and F_3 separately, without floating point.
These finite tests do not prove the universal assertion or field independence.
"""
from collections import deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product

FIELDS = (0, 2, 3)  # 0 denotes Q.


def graph(n, mask):
    adj = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if mask & (1 << bit):
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    return tuple(adj)


def subgraph(adj, subset):
    vs = [v for v in range(len(adj)) if subset & (1 << v)]
    return tuple(sum(1 << j for j, v in enumerate(vs) if adj[u] & (1 << v))
                 for u in vs)


def independent(adj, subset):
    return all(not (adj[v] & subset) for v in range(len(adj)) if subset >> v & 1)


@lru_cache(None)
def power(adj, radius):
    out = []
    for start in range(len(adj)):
        distances = {start: 0}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in range(len(adj)):
                if adj[u] >> v & 1 and v not in distances:
                    distances[v] = distances[u] + 1
                    queue.append(v)
        out.append(sum(1 << v for v, d in distances.items() if 0 < d <= radius))
    return tuple(out)


def rank(rows, characteristic):
    if not rows:
        return 0
    a = [[Fraction(x) if characteristic == 0 else x % characteristic for x in row]
         for row in rows]
    lead = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(lead, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[lead], a[pivot] = a[pivot], a[lead]
        inv = 1 / a[lead][col] if characteristic == 0 else pow(a[lead][col], -1, characteristic)
        a[lead] = [x * inv for x in a[lead]]
        if characteristic:
            a[lead] = [x % characteristic for x in a[lead]]
        for i in range(lead + 1, len(a)):
            coef = a[i][col]
            if coef:
                a[i] = [x - coef * y for x, y in zip(a[i], a[lead])]
                if characteristic:
                    a[i] = [x % characteristic for x in a[i]]
        lead += 1
        if lead == len(a):
            break
    return lead


@lru_cache(None)
def reduced_homology(adj, characteristic):
    n = len(adj)
    faces = [[s for s in range(1 << n) if s.bit_count() == r and independent(adj, s)]
             for r in range(n + 1)]
    ranks = [0] * (n + 2)
    for r in range(1, n + 1):
        lower_index = {s: i for i, s in enumerate(faces[r - 1])}
        matrix = [[0] * len(faces[r]) for _ in faces[r - 1]]
        for j, s in enumerate(faces[r]):
            vertices = [v for v in range(n) if s >> v & 1]
            for i, v in enumerate(vertices):
                matrix[lower_index[s ^ (1 << v)]][j] = (-1) ** i
        ranks[r] = rank(matrix, characteristic)
    betti = tuple(len(faces[r]) - ranks[r] - ranks[r + 1] for r in range(n + 1))
    assert all(b >= 0 for b in betti)
    return betti  # entry r = dim reduced H_(r-1).


@lru_cache(None)
def pd_quotient(adj, characteristic):
    value = 0
    for subset in range(1 << len(adj)):
        w = subset.bit_count()
        for r, h in enumerate(reduced_homology(subgraph(adj, subset), characteristic)):
            if h:
                value = max(value, w - r)
    return value


def pd_ideal(adj, characteristic):
    if not any(adj):
        return None  # bottom / negative infinity, as in Mathlib.
    result = pd_quotient(adj, characteristic) - 1
    assert result >= 0
    return result


def alpha(adj):
    return max(s.bit_count() for s in range(1 << len(adj)) if independent(adj, s))


def tau(adj):
    assert all(adj)  # no isolated vertices: every domination minimum is nonempty.
    n = len(adj)
    unions = [0] * (1 << n)
    for x in range(1, 1 << n):
        bit = x & -x
        unions[x] = unions[x ^ bit] | adj[bit.bit_length() - 1]
    return max(min(x.bit_count() for x in range(1 << n) if not (a & ~unions[x]))
               for a in range(1 << n) if independent(adj, a))


def verify(adj):
    n = len(adj)
    square = power(adj, 2)
    for characteristic in FIELDS:
        lhs = pd_ideal(adj, characteristic)
        for r in range(2, max(3, n)):
            rhs = pd_ideal(power(adj, r), characteristic)
            assert (lhs is None and rhs is None) or (lhs is not None and rhs is not None and lhs <= rhs), (adj, r, characteristic, lhs, rhs)
        if n and all(adj):
            t = tau(adj)
            assert t >= alpha(square)
            assert pd_quotient(adj, characteristic) <= n - t <= n - alpha(square) <= pd_quotient(square, characteristic)


def prufer_tree(n, word):
    degrees = [1] * n
    for x in word:
        degrees[x] += 1
    adj = [0] * n
    for x in word:
        v = next(v for v in range(n) if degrees[v] == 1)
        adj[x] |= 1 << v
        adj[v] |= 1 << x
        degrees[x] -= 1
        degrees[v] -= 1
    if n >= 2:
        u, v = [v for v in range(n) if degrees[v] == 1]
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return tuple(adj)


def main():
    fidelity = []
    for n in range(2, 9):
        adj = tuple(sum(1 << v for v in range(n) if abs(u - v) == 1) for u in range(n))
        expected_square = n - (n + 4) // 5
        expected_path = n - (n + 2) // 3
        for characteristic in FIELDS:
            assert pd_quotient(power(adj, 2), characteristic) == expected_square
            assert pd_quotient(adj, characteristic) == expected_path
        fidelity.append((n, expected_path, expected_square))
    print('FIDELITY: source Proposition 4.3 and Corollary 4.4 (v2), quotient pd(P_n^2)=n-ceil(n/5), 2<=n<=8')
    print('PATH_ROWS (n, pd S/I(P_n), pd S/I(P_n^2)):', fidelity)
    graph_count = 0
    for n in range(6):
        count = 1 << (n * (n - 1) // 2)
        for mask in range(count):
            verify(graph(n, mask))
        graph_count += count
        print('SIMPLE_GRAPHS n=%d count=%d OK' % (n, count))
    tree_count = 0
    for n in range(1, 7):
        count = 0
        for word in product(range(n), repeat=max(0, n - 2)):
            verify(prufer_tree(n, word))
            count += 1
        tree_count += count
        print('LABELLED_TREES n=%d count=%d OK' % (n, count))
    print('CHECKED: %d labelled simple graphs (0<=n<=5); %d labelled trees (1<=n<=6)' % (graph_count, tree_count))
    print('FIELDS: Q, F_2, F_3; distance powers 2<=r<=max(2,n-1); zero ideals and disconnected graphs included')
    print('METHOD: Hochster formula; exact reduced boundary ranks; no floating point')
    print('LIMIT: finite evidence only; neither a universal proof nor a check over every field')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
