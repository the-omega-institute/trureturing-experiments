#!/usr/bin/env python3
"""Exact checks for balanced complete multipartite general-position counts.

Source: B. A. Rather, arXiv:2603.06930v3, Sections 3--4.
The proposed classification is checked on a finite parameter rectangle.
These computations are corroboration of a proposed uniform argument.
"""
from collections import deque
from fractions import Fraction
from itertools import combinations
from math import comb


def binomial_row(n):
    row = [1]
    for k in range(n):
        row.append(row[-1] * (n - k) // (k + 1))
    return row


def coefficients(t, n):
    at, bn = binomial_row(t), binomial_row(n)
    return [1, t * n] + [
        (n**k * at[k] if k <= t else 0) + (t * bn[k] if k <= n else 0)
        for k in range(2, max(t, n) + 1)
    ]


def log_concave(a):
    return all(a[k] ** 2 >= a[k - 1] * a[k + 1] for k in range(1, len(a) - 1))


def unimodal(a):
    # Independent implementation of weak increase followed by weak decrease.
    falling = False
    for left, right in zip(a, a[1:]):
        if right < left:
            falling = True
        if falling and right > left:
            return False
    return True


def predicted_lc(t, n):
    return n <= t + 1 or (t, n) == (2, 4)


def junction_ratio(t, n):
    return Fraction(t * comb(n, t) * (n - 2 * t - 1), (t + 1) * n**t)


def predicted_unimodal(t, n):
    return n <= 2 * t + 3 or junction_ratio(t, n) >= 1


def graph_counts(t, n):
    """Enumerate subsets using BFS graph distances, without the counting formula."""
    order = t * n
    adjacency = [[v for v in range(order) if u // n != v // n] for u in range(order)]
    distance = []
    for start in range(order):
        dist = [-1] * order
        dist[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        assert all(d >= 0 for d in dist)
        distance.append(dist)
    forbidden = []
    for u, v, w in combinations(range(order), 3):
        if (distance[u][w] == distance[u][v] + distance[v][w]
                or distance[u][v] == distance[u][w] + distance[w][v]
                or distance[v][w] == distance[v][u] + distance[u][w]):
            forbidden.append((1 << u) | (1 << v) | (1 << w))
    counts = [0] * (order + 1)
    for mask in range(1 << order):
        if all(mask & triple != triple for triple in forbidden):
            counts[mask.bit_count()] += 1
    while counts[-1] == 0:
        counts.pop()
    return counts


def check_fidelity():
    # Source Section 4 examples, including a strict valley and LC-only failure.
    assert coefficients(2, 8) == [1, 16, 120, 112, 140, 112, 56, 16, 2]
    assert coefficients(2, 5) == [1, 10, 45, 20, 10, 2]
    assert not unimodal(coefficients(2, 8))
    assert unimodal(coefficients(2, 5)) and not log_concave(coefficients(2, 5))
    for t in range(2, 101):
        for n in range(1, 5):
            a = coefficients(t, n)
            assert log_concave(a) and unimodal(a)
    print('FIDELITY: source K_(8,8), K_(5,5) coefficient lists reproduced')
    print('FIDELITY: source small-part theorem checked for 2 <= t <= 100, 1 <= n <= 4 (396 pairs)')


def check_graphs():
    cases = [(t, n) for t in range(2, 5) for n in range(1, 5)] + [(5, 2), (5, 3)]
    total_subsets = 0
    for t, n in cases:
        a = graph_counts(t, n)
        assert a == coefficients(t, n), (t, n, a)
        # The polynomial identity, including the empty-set/singleton correction.
        for k in range(t * n + 1):
            rhs = (comb(t, k) * n**k if k <= t else 0)
            rhs += t * (comb(n, k) if k <= n else 0)
            rhs -= t if k == 0 else t * n if k == 1 else 0
            assert rhs == (a[k] if k < len(a) else 0)
        total_subsets += 1 << (t * n)
    print(f'GRAPH_COUNTS: BFS geodesics and all subsets agree on {len(cases)} graphs; {total_subsets} subsets; largest order 16')


def check_classification():
    lc_pairs = um_pairs = non_um_pairs = 0
    for t in range(2, 16):
        for n in range(1, 401):
            a = coefficients(t, n)
            lc, um = log_concave(a), unimodal(a)
            assert lc == predicted_lc(t, n), ('LC', t, n)
            assert um == predicted_unimodal(t, n), ('UM', t, n)
            failures = [k for k in range(1, len(a) - 1) if a[k] ** 2 < a[k-1] * a[k+1]]
            assert all(k == t + 1 for k in failures), ('LC failure location', t, n, failures)
            if n >= t + 2:
                assert all(a[k] < a[k+1] for k in range(t)), ('prefix', t, n)
            if n >= 2 * t + 4:
                assert um == (a[t] <= a[t+1]), ('junction', t, n)
            lc_pairs += lc
            um_pairs += um
            non_um_pairs += not um
    print('CLASSIFICATION: 2 <= t <= 15, 1 <= n <= 400: 5600 pairs agree')
    print(f'CLASSIFICATION: log-concave={lc_pairs}; unimodal={um_pairs}; nonunimodal={non_um_pairs}')
    print('STRUCTURE: every checked LC failure is at t+1; increasing prefix and junction test agree')


def check_thresholds():
    thresholds = []
    for t in range(2, 16):
        start = 2 * t + 4
        assert junction_ratio(t, start) < 1
        a = coefficients(t, start)
        assert a[t] > a[t+1] < a[t+2], ('family valley', t)
        lo, hi = start, 2 * start
        while junction_ratio(t, hi) < 1:
            hi *= 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if junction_ratio(t, mid) >= 1:
                hi = mid
            else:
                lo = mid
        assert junction_ratio(t, hi - 1) < 1 <= junction_ratio(t, hi)
        thresholds.append((t, hi))
    print('THRESHOLD_BOUNDARIES: exact junction tests; monotonicity is a separate proof obligation')
    print('N_t=' + ', '.join(f'{t}:{nt}' for t, nt in thresholds))
    print('FAMILY: n=2t+4 has a strict valley at t+1 for 2 <= t <= 15')


def main():
    check_fidelity()
    check_graphs()
    check_classification()
    check_thresholds()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
