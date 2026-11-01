#!/usr/bin/env python3
"""Exact finite checks for arXiv:2608.10491v1, Conjecture 1.
The companion shift follows the source's concatenated history definition.
Finite experimental evidence is not a proof for arbitrary parameters.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from collections import Counter


def solve(a, b):
    """Gauss-Jordan over Q, with exact independent residual checks."""
    n = len(b)
    rows = [[Q(x) for x in row] + [Q(y)] for row, y in zip(a, b)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if rows[i][k])
        rows[k], rows[pivot] = rows[pivot], rows[k]
        d = rows[k][k]
        rows[k] = [x / d for x in rows[k]]
        for i in range(n):
            if i != k:
                d = rows[i][k]
                rows[i] = [x - d * y for x, y in zip(rows[i], rows[k])]
    x = tuple(row[-1] for row in rows)
    assert all(sum(u * v for u, v in zip(row, x)) == y for row, y in zip(a, b))
    return x


def matrices(n):
    rows = [r for r in product((Q(0), Q(1, 2), Q(1)), repeat=n) if sum(r) == 1]
    return product(rows, repeat=n)


def equilibrium(s, c, d, delays):
    n = len(c)
    depth = max(max(x) for x in delays)
    size = n * depth
    j = [[Q(0) for _ in range(size)] for _ in range(size)]
    for i in range(n):
        for k in range(n):
            for lag in delays[i * n + k]:
                j[i][(lag - 1) * n + k] += d[i] * s[i][k]
    for lag in range(1, depth):
        for i in range(n):
            j[lag * n + i][(lag - 1) * n + i] = Q(1)
    b = tuple(c) + (Q(0),) * (n * (depth - 1))
    a = [[Q(i == k) - j[i][k] for k in range(size)] for i in range(size)]
    x = solve(a, b)
    p = x[:n]
    assert all(x[lag * n + i] == p[i] for lag in range(depth) for i in range(n))
    return p


def reduced(s, c, d, m):
    n = len(c)
    a = [[Q(i == k) - m * d[i] * s[i][k] for k in range(n)] for i in range(n)]
    return solve(a, c)


def delay_arrays(n, t, m):
    return product(tuple(combinations(range(1, t + 1), m)), repeat=n * n)


def fidelity():
    # Source eqn:homo_fixed_point: p* = c_{1,M}/(1-M c_{2,M}).
    counts = Counter()
    for n in (1, 2):
        for s in matrices(n):
            for t in (1, 2, 3):
                for m in range(1, t + 1):
                    c, d = (Q(1, 4),) * n, (Q(1, 8),) * n
                    expected = Q(1, 4) / (1 - m * Q(1, 8))
                    for delays in delay_arrays(n, t, m):
                        assert equilibrium(s, c, d, delays) == (expected,) * n
                        counts[n] += 1
    print('FIDELITY homogeneous closed form:', dict(sorted(counts.items())))


def exhaustive_small():
    counts = Counter()
    minimum_gain = None
    for n in (1, 2):
        for s in matrices(n):
            for c in product((Q(1, 8), Q(1, 4)), repeat=n):
                for d in product((Q(1, 16), Q(1, 8)), repeat=n):
                    for t in (1, 2, 3):
                        solutions = {}
                        for m in range(1, t + 1):
                            assert all(0 < c[i] and 0 < d[i] and c[i] + m * d[i] < 1
                                       for i in range(n))
                            expected = reduced(s, c, d, m)
                            for delays in delay_arrays(n, t, m):
                                p = equilibrium(s, c, d, delays)
                                assert p == expected
                                assert all(0 < z < 1 for z in p)
                                counts['augmented_systems'] += 1
                            solutions[m] = expected
                        for m1, m2 in combinations(range(1, t + 1), 2):
                            gains = tuple(y - x for x, y in zip(solutions[m1], solutions[m2]))
                            assert all(z > 0 for z in gains)
                            minimum_gain = min(gains) if minimum_gain is None else min(minimum_gain, *gains)
                            counts['memory_pairs'] += 1
    print('SMALL DOMAIN: N=1,2; T=1,2,3; all stochastic half-grid rows;')
    print('  c_i in {1/8,1/4}; d_i in {1/16,1/8}; every per-link M-subset of {1,...,T}.')
    print('COUNTS:', dict(sorted(counts.items())))
    print('MINIMUM STRICT COORDINATE GAIN:', minimum_gain)
    print('DELAY COMPARISON: every architecture at each M has the same exact solution;')
    print('  hence every cross-architecture pair in this finite domain is compared.')


def larger_examples():
    n, t = 3, 5
    networks = (
        ((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1))),
        ((Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1)), (Q(1), Q(0), Q(0))),
        ((Q(1, 2), Q(1, 2), Q(0)), (Q(0), Q(1), Q(0)), (Q(1, 3), Q(1, 3), Q(1, 3))),
    )
    c, d = (Q(1, 10), Q(1, 5), Q(3, 10)), (Q(1, 20), Q(1, 16), Q(1, 12))
    count = 0
    for s in networks:
        old = None
        for m in range(1, t + 1):
            choices = tuple(combinations(range(1, t + 1), m))
            expected = reduced(s, c, d, m)
            for seed in range(7):
                delays = tuple(choices[(seed + 3 * i + 7 * k) % len(choices)]
                               for i in range(n) for k in range(n))
                assert equilibrium(s, c, d, delays) == expected
                count += 1
            assert all(0 < x < 1 for x in expected)
            if old is not None:
                assert all(x < y for x, y in zip(old, expected))
            old = expected
    print('LARGER EXAMPLES: N=3, T=5; identity, cycle, reducible networks;', count, 'exact systems.')


if __name__ == '__main__':
    fidelity()
    exhaustive_small()
    larger_examples()
    print('ALL_CHECKS_OK')
