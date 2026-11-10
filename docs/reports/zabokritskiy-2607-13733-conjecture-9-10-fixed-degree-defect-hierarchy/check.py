#!/usr/bin/env python3
"""Exact finite checks for arXiv:2607.13733v3, Conjecture 9.10.

No floating point, external packages, random inputs or network access.
These finite checks do not prove the unbounded conjecture.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from itertools import permutations
from math import comb, factorial, prod


@cache
def partitions(n, largest=None):
    if largest is None:
        largest = n
    if n == 0:
        return ((),)
    return tuple((a,) + tail for a in range(min(n, largest), 0, -1)
                 for tail in partitions(n - a, a))


def defect(lam):
    return sum(lam) - len(lam)


def predecessor(mu, s):
    assert s in mu and s > 0
    a = list(mu)
    a.remove(s)
    if s > 1:
        a.append(s - 1)
    return tuple(sorted(a, reverse=True))


@cache
def derangements(n):
    # Inclusion-exclusion; checked independently by direct permutations below.
    return sum((-1)**j * factorial(n) // factorial(j) for j in range(n + 1))


@cache
def fixed_counts(n):
    return tuple(comb(n, j) * derangements(n - j) for j in range(n + 1))


@cache
def total_fixed_counts(lam):
    counts = [1]
    for a in lam:
        block = fixed_counts(a)
        nxt = [0] * (len(counts) + a)
        for f, u in enumerate(counts):
            for g, v in enumerate(block):
                nxt[f + g] += u * v
        counts = nxt
    assert sum(counts) == prod(factorial(a) for a in lam)
    return tuple(counts)


@cache
def average_weight(lam):
    counts = total_fixed_counts(lam)
    return sum((Q(c, 1 + f) for f, c in enumerate(counts)), Q(0)) / sum(counts)


def recursive_weights(d):
    # The source's exact-coverage equation at mu = lam union {1}:
    # (m_1(lam)+1) w_lam + sum_{s>=2} m_s(lam) w_pred_s(mu) = 1.
    # Each weight in the sum has strictly smaller defect.
    w = {}
    for lam in sorted(partitions(d - 1), key=defect):
        mu = tuple(sorted(lam + (1,), reverse=True))
        multiplicities = Counter(mu)
        lower = sum((m * w[predecessor(mu, s)]
                     for s, m in multiplicities.items() if s > 1), Q(0))
        w[lam] = (1 - lower) / multiplicities[1]
    return w


def coverage(mu, w):
    return sum((m * w[predecessor(mu, s)] for s, m in Counter(mu).items()), Q(0))


def main():
    direct = 0
    for n in range(8):
        observed = Counter(sum(i == v for i, v in enumerate(p))
                           for p in permutations(range(n)))
        assert tuple(observed[j] for j in range(n + 1)) == fixed_counts(n)
        direct += factorial(n)
    print(f'Fixed-point distributions: direct enumeration of {direct} permutations, 0<=n<=7: OK')

    source_patterns = ((1,) * 6, (2,) + (1,) * 4,
                       (3,) + (1,) * 3, (2, 2) + (1,) * 2)
    expected = (Q(1, 7), Q(6, 35), Q(29, 140), Q(23, 105))
    rec7 = recursive_weights(7)
    assert tuple(rec7[p] for p in source_patterns) == expected
    assert tuple(average_weight(p) for p in source_patterns) == expected
    print('Fidelity: Eq. (69), d=7: ' + ', '.join(map(str, expected)) + ': OK')

    compared = 0
    for d in range(1, 15):
        for lam, w in recursive_weights(d).items():
            assert w == average_weight(lam), (d, lam, w, average_weight(lam))
            assert 0 < w <= 1
            compared += 1
    print(f'Recursion versus independent fixed-point average: {compared} anchor partitions, 1<=d<=14: OK')

    identities = 0
    for d in range(1, 17):
        for mu in partitions(d):
            lhs = sum((m * average_weight(predecessor(mu, s))
                       for s, m in Counter(mu).items()), Q(0))
            counts = total_fixed_counts(mu)
            rhs = 1 - Q(counts[0], sum(counts))
            assert lhs == rhs, (mu, lhs, rhs)
            if 1 in mu:
                assert lhs == 1
            identities += 1
    print(f'Marked-fixed-point identity: {identities} vertex partitions, 1<=d<=16: OK')

    checked = low = nonnegative = levels = 0
    for d in range(3, 17):
        recursive = recursive_weights(d)
        for k in range((d - 3) // 2 + 1):
            assert 2 * k + 3 <= d
            levels += 1
            w = {lam: value if defect(lam) <= k else Q(1)
                 for lam, value in recursive.items()}
            for lam, value in w.items():
                assert value >= 0, (d, k, lam, value)
                nonnegative += 1
            for mu in partitions(d):
                c = coverage(mu, w)
                assert c >= 1, (d, k, mu, c)
                checked += 1
                if defect(mu) <= k:
                    assert c == 1, (d, k, mu, c)
                    low += 1
    print(f'Truncated hierarchy: {levels} (d,k) pairs, 2k+3<=d<=16: OK')
    print(f'Nonnegative anchor weights: {nonnegative}; coverage inequalities: {checked}; exact low-defect coverages: {low}: OK')

    for k in range(8):
        d = 2 * k + 2
        w = {lam: v if defect(lam) <= k else Q(1)
             for lam, v in recursive_weights(d).items()}
        mu = (2,) * (k + 1)
        assert coverage(mu, w) == 1 - Q(1, 2**(k + 1)) < 1
    print('Triage only: threshold d=2k+2 fails on (2^(k+1)), 0<=k<=7: OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
