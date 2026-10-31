#!/usr/bin/env python3
"""Exact finite checks for Kurkov's first A140717 FORMULA conjecture.

Only integers and tuples of integer polynomial coefficients are used.
Finite checks support a preregistration; they are not a proof for all sizes.
"""
from collections import Counter
from itertools import combinations, permutations
from math import comb

ZERO = (0,)
ONE = (1,)


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    while len(c) > 1 and not c[-1]:
        c.pop()
    return tuple(c)


def shift(a, k):
    assert k >= 0
    return (0,) * k + a if a != ZERO else ZERO


def histogram(weights):
    h = Counter(weights)
    return tuple(h[k] for k in range(max(h, default=0) + 1))


def triangular(m):
    assert m >= 1
    v = [ONE] * m
    t = v.copy()
    stages = {}
    for i in range(m - 1):
        for j in range(i + 1, m):
            a, b = v[i], v[j]
            v[i], v[j] = add(a, shift(b, j - i)), add(shift(a, 1), b)
        t[i + 1] = v[i + 1]
        for j in range(i + 1, m):
            stages[i + 1, j + 1] = v[j]
    return t, stages


def dyck_words(n, up=0, down=0, prefix=''):
    if down == n:
        yield prefix
    else:
        if up < n:
            yield from dyck_words(n, up + 1, down, prefix + 'U')
        if down < up:
            yield from dyck_words(n, up, down + 1, prefix + 'D')


def peak_statistic(word):
    height = 0
    weight = 0
    for i, step in enumerate(word):
        height += 1 if step == 'U' else -1
        assert height >= 0
        if step == 'U' and i + 1 < len(word) and word[i + 1] == 'D':
            weight += height - 1
    assert height == 0
    return weight


def peak_pairs(word):
    """A peak after u up-steps and d down-steps gives (u-1,d+1)."""
    up = down = 0
    result = []
    for i, step in enumerate(word):
        up += step == 'U'
        down += step == 'D'
        if step == 'U' and i + 1 < len(word) and word[i + 1] == 'D' and up - down >= 2:
            result.append((up - 1, down + 1))
    return tuple(result)


def pair_weight(pairs):
    return sum(a - b + 1 for a, b in pairs)


def pair_sequences(i, j):
    """Strictly increasing a<j and b<=i, with 1<=b<=a."""
    for r in range(min(i, j - 1) + 1):
        for aa in combinations(range(1, j), r):
            for bb in combinations(range(1, i + 1), r):
                if all(b <= a for a, b in zip(aa, bb)):
                    yield tuple(zip(aa, bb))


def decode(n, pairs):
    """Reconstruct a Dyck path from its nontrivial peaks, filling at baseline."""
    word = ''
    x = y = 0
    for a, b in pairs:
        u, d = a + 1, b - 1
        assert y <= d and x < u and d < u <= n
        if x <= d:
            word += 'D' * (x - y)
            word += 'UD' * (d - x)
            word += 'U' * (u - d)
        else:
            word += 'D' * (d - y)
            word += 'U' * (u - x)
        x, y = u, d
    word += 'D' * (x - y)
    word += 'UD' * (n - x)
    return word


def rectangle_recurrence(max_n):
    f = {(0, j): ONE for j in range(1, max_n + 1)}
    for i in range(1, max_n):
        for j in range(i + 1, max_n + 1):
            value = f[i - 1, j]
            for a in range(i, j):
                value = add(value, shift(f[i - 1, a], a - i + 1))
            f[i, j] = value
    return f


def avoids321(p):
    return not any(p[i] > p[j] > p[k]
                   for i in range(len(p)) for j in range(i + 1, len(p))
                   for k in range(j + 1, len(p)))


def inversions(p):
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def main():
    # OEIS A140717 EXAMPLE rows, starting at semilength zero.
    source_rows = [(1,), (1,), (1, 1), (1, 2, 2), (1, 3, 5, 4, 1),
                   (1, 4, 9, 12, 10, 4, 2), (1, 5, 14, 25, 31, 26, 16, 9, 4, 1)]
    dyck_rows = {}
    total_paths = 0
    for n in range(11):
        words = list(dyck_words(n))
        assert len(words) == comb(2 * n, n) // (n + 1)
        total_paths += len(words)
        dyck_rows[n] = histogram(peak_statistic(w) for w in words)
        assert len(dyck_rows[n]) - 1 == n * n // 4
        for w in words:
            pairs = peak_pairs(w)
            assert all(1 <= b <= a < n for a, b in pairs)
            assert all(a < c and b < d for (a, b), (c, d) in zip(pairs, pairs[1:]))
            assert pair_weight(pairs) == peak_statistic(w)
            assert decode(n, pairs) == w
        for pairs in pair_sequences(max(n - 1, 0), n or 1):
            w = decode(n, pairs)
            assert peak_pairs(w) == pairs
            assert peak_statistic(w) == pair_weight(pairs)
        if n < len(source_rows):
            assert dyck_rows[n] == source_rows[n]
    print('FIDELITY: OEIS printed rows n=0..6 reproduced exactly.')
    print('FIDELITY: T(4,2)=5; R(4,z) coefficients =', dyck_rows[4])
    print('DYCK_BIJECTION: both inverse identities and weight preservation for n=0..10;', total_paths, 'paths.')
    f = rectangle_recurrence(30)
    for j in range(1, 11):
        for i in range(j):
            enumerated = histogram(pair_weight(p) for p in pair_sequences(i, j))
            assert enumerated == f[i, j]
    print('RECTANGLE_RECURRENCE: direct pair enumeration agrees for 0<=i<j<=10.')
    largest, _ = triangular(30)
    for m in range(1, 31):
        outputs, stages = triangular(m)
        assert outputs == largest[:m]
        for ij, value in stages.items():
            assert value == f[ij]
        for n in range(1, min(m, 10) + 1):
            assert outputs[n - 1] == dyck_rows[n]
        for n in range(1, m + 1):
            assert sum(outputs[n - 1]) == comb(2 * n, n) // (n + 1)
            assert len(outputs[n - 1]) - 1 == n * n // 4
    print('ALGORITHM: all saved outputs and intermediate states agree with the recurrence for every m=1..30.')
    print('CLAIM_SMALL_CASES: every 1<=n<=m<=30 with n<=10 agrees with direct Dyck enumeration; all prefixes stable.')
    for n in range(9):
        row = histogram(inversions(p) for p in permutations(range(1, n + 1)) if avoids321(p))
        assert row == dyck_rows[n]
    print('INDEPENDENT_CROSSCHECK: 321-avoiding permutations, n=0..8, agree by inversion number.')
    # Two source examples are erroneous; neither is part of the target formula.
    actual = [p for p in permutations(range(1, 5)) if avoids321(p) and inversions(p) == 2]
    assert len(actual) == 5 and inversions((2, 3, 4, 1)) == 3
    assert len('UDUUUDD') == 7
    print('SOURCE_ERRATA: printed permutation 2341 has 3 inversions; replace it by 2314 in the T(4,2) list.')
    print('SOURCE_ERRATA: printed word UDUUUDD has length 7; UDUUUDDD is a valid missing semilength-4 example.')
    assert peak_statistic('UDUUUDDD') == 2
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
