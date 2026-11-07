#!/usr/bin/env python3
"""Exact checks for Conjecture 1 of arXiv:1812.09458v2. Python standard library only."""
from fractions import Fraction
from itertools import combinations, product
from math import comb


def mul(a, b):
    n = len(a)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for j, x in enumerate(a[i]):
            if x:
                for k, y in enumerate(b[j]):
                    out[i][k] += x * y
    return out


def laplacian(a):
    return [[sum(a[i]) if i == j else -a[i][j]
             for j in range(len(a))] for i in range(len(a))]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def norm(a):
    return max((sum(abs(x) for x in row) for row in a), default=0)


def tournament(n, mask):
    a = [[0] * n for _ in range(n)]
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if mask & (1 << bit):
            a[i][j] = 1
        else:
            a[j][i] = 1
    return a


def check_tournament(a):
    n = len(a)
    assert all(a[i][i] == 0 for i in range(n))
    assert all(a[i][j] in (0, 1) and a[i][j] + a[j][i] == 1
               for i, j in combinations(range(n), 2))
    assert sum(map(sum, a)) == comb(n, 2)


def closed_walk_trace(l, k):
    """A separate trace computation by cyclic products of entries."""
    total = 0
    for vertices in product(range(len(l)), repeat=k):
        value = 1
        for i in range(k):
            value *= l[vertices[i]][vertices[(i + 1) % k]]
        total += value
    return total


def source_h2_formula(n):
    if n % 2:
        return Fraction(comb(n + 1, 3), 4) + 1
    return 2 * comb(n // 2 + 1, 3) + 1


def enumerate_small():
    expected_s = {2: 1, 3: 2, 4: 4, 5: 9, 6: 22}
    expected_h2 = {2: 1, 3: 2, 4: 3, 5: 6, 6: 9}
    for n in range(2, 7):
        scores = set()
        readings = {k: set() for k in range(2, 9)}
        positive = {k: set() for k in range(2, 9)}
        for mask in range(1 << comb(n, 2)):
            a = tournament(n, mask)
            check_tournament(a)
            degrees = tuple(sorted(map(sum, a)))
            scores.add(degrees)
            l = laplacian(a)
            power = l
            for k in range(2, 9):
                power = mul(power, l)
                t = trace(power)
                bound = n * (2 * (n - 1)) ** k
                assert norm(power) <= (2 * (n - 1)) ** k
                assert abs(t) <= bound
                assert Fraction(t, comb(n, 2) ** k) == trace(
                    [[Fraction(x, comb(n, 2) ** k) for x in row] for row in power])
                if k == 2:
                    assert t == sum(d * d for d in degrees)
                    assert t > 0
                if n <= 4 and k <= 4:
                    assert t == closed_walk_trace(l, k)
                readings[k].add(t)
                if t > 0:
                    positive[k].add(t)
        assert len(scores) == expected_s[n]
        assert len(readings[2]) == expected_h2[n] == source_h2_formula(n)
        if n == 6:
            assert len(readings[4]) == 36
        print('EXHAUSTIVE', 'n=' + str(n), 'labeled=' + str(1 << comb(n, 2)),
              'S=' + str(len(scores)), 'h2=' + str(len(readings[2])),
              'trace_counts_k2_to_k8=' + str([len(readings[k]) for k in readings]),
              'positive_trace_counts=' + str([len(positive[k]) for k in positive]))


def block_family(bits, residue):
    """Later blocks beat earlier blocks; bit 0 is transitive, bit 1 is a triangle."""
    m = len(bits)
    n = 3 * m + residue
    block_of = [i // 3 if i < 3 * m else m + i - 3 * m for i in range(n)]
    a = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if block_of[i] != block_of[j]:
            a[j][i] = 1
        elif bits[block_of[i]]:
            if (j - i) % 3 == 1:
                a[i][j] = 1
            else:
                a[j][i] = 1
        else:
            a[j][i] = 1
    return a


def check_family():
    total = 0
    for m in range(9):
        for residue in range(3):
            distinct = set()
            for bits in product((0, 1), repeat=m):
                a = block_family(bits, residue)
                check_tournament(a)
                scores = tuple(sorted(map(sum, a)))
                expected = []
                for j, bit in enumerate(bits):
                    expected.extend([3*j + 1] * 3 if bit else [3*j, 3*j + 1, 3*j + 2])
                expected.extend(range(3*m, 3*m + residue))
                assert scores == tuple(expected)
                recovered = tuple(int(scores[3*j:3*j+3] == (3*j+1,) * 3) for j in range(m))
                assert recovered == bits
                distinct.add(scores)
                total += 1
            assert len(distinct) == 2 ** m
    print('SCORE_FAMILY m=0..8 residues=0,1,2 realized=' + str(total),
          'distinct_for_each_order=2^m; every bit recovered from actual outdegrees')


def check_positivity_bound():
    for k in range(2, 13):
        n = 8 ** k + 1
        lower = n ** k * (Fraction(n, 4 ** k) - (2 ** k - 1))
        assert lower > 0
        assert Fraction(n - 1, 2) >= Fraction(n, 4)
    print('POSITIVITY_LOWER_BOUND k=2..12 n=8^k+1 exact rational lower bound > 0')


def check_tail_majorant():
    n = 900
    for k in range(2, 13):
        bound = 2 * n * (2 * (n - 1)) ** k + 1
        upper = Fraction(bound, 2 ** (n // 3))
        assert upper * 2 ** 64 < 1
    print('TAIL_MAJORANT n=900 k=2..12: (2*n*(2*(n-1))^k+1)/2^floor(n/3) < 1/2^64')


def check_domain_edge():
    cycle = block_family((1,), 0)
    p = laplacian(cycle)
    values = {}
    for k in range(2, 7):
        p = mul(p, laplacian(cycle))
        values[k] = trace(p)
    assert values[3] == 0 and values[4] < 0
    print('DOMAIN_EDGE C3 raw traces k=2..6=' + str([values[k] for k in values]),
          '; nonpositive small-order moments are excluded from real entropy values')


if __name__ == '__main__':
    enumerate_small()
    check_family()
    check_positivity_bound()
    check_tail_majorant()
    check_domain_edge()
    print('ALL_CHECKS_OK')
