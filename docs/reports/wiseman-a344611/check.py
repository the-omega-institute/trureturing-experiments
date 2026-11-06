#!/usr/bin/env python3
"""Exact finite checks for the A344611 preregistration; no floating-point evidence.

Sources: https://oeis.org/A344611 (revision 27), A160786 (revision 25),
A027193 (revision 78), A344743 (revision 10), A236914 (revision 23).
Enumeration is independently cross-checked by an unbounded-parts parity DP.
The finite checks do not prove the universally quantified conjecture.
"""
from collections import Counter

MAX_N = 20
SOURCE_A = (1, 2, 4, 8, 15, 27, 48, 81, 135, 220, 352, 553, 859, 1313, 1986, 2969)
SOURCE_B = (1, 2, 4, 8, 16, 29, 52, 90, 151, 248, 400, 632, 985, 1512, 2291, 3431)
SOURCE_C = (0, 0, 1, 3, 7, 15, 29, 54, 96, 165, 275, 449, 716, 1123, 1732, 2635)
SOURCE_D = (0, 1, 3, 7, 14, 27, 49, 86, 146, 242, 392, 623, 973, 1498, 2274, 3411)
SOURCE_DIFF = (0, 0, 0, 0, 1, 2, 4, 9, 16, 28, 48, 79)
SOURCE_EXAMPLE_WEIGHT_8 = {
    (8,), (4, 4), (3, 3, 2), (4, 2, 2), (4, 3, 1), (5, 2, 1), (6, 1, 1),
    (2, 2, 2, 2), (3, 3, 1, 1), (2, 2, 2, 1, 1), (3, 2, 1, 1, 1),
    (4, 1, 1, 1, 1), (2, 2, 1, 1, 1, 1), (2, 1, 1, 1, 1, 1, 1), (1,) * 8,
}


def partitions(weight, ceiling=None):
    """Every nonincreasing positive tuple of the specified weight, once."""
    if weight == 0:
        yield ()
        return
    ceiling = weight if ceiling is None else min(weight, ceiling)
    for first in range(ceiling, 0, -1):
        for tail in partitions(weight - first, first):
            yield (first,) + tail


def reverse_sum(parts):
    """The source formula, i=1..k: sum (-1)^(k-i) y_i."""
    k = len(parts)
    return sum((-1 if (k - i - 1) % 2 else 1) * value
               for i, value in enumerate(parts))


def recurrence_sum(parts):
    endpoint = 0
    for value in parts:
        endpoint = value - endpoint
    return endpoint


def conjugate(parts):
    return tuple(sum(value >= column for value in parts)
                 for column in range(1, parts[0] + 1)) if parts else ()


def injection(parts):
    if len(parts) % 2:
        return (parts[0] + 1,) + parts[1:]
    return parts + (1,)


def parity_dp(bound):
    """Independent count by weight and parity of length, using unbounded parts."""
    counts = [[0, 0] for _ in range(bound + 1)]
    counts[0][0] = 1
    for part in range(1, bound + 1):
        for weight in range(part, bound + 1):
            even, odd = counts[weight - part]
            counts[weight][0] += odd
            counts[weight][1] += even
    return counts


def main():
    dp = parity_dp(2 * MAX_N + 1)
    rows = []
    checked_sources = 0
    for n in range(MAX_N + 1):
        even_weight = list(partitions(2 * n))
        odd_weight = list(partitions(2 * n + 1))
        assert len(set(even_weight)) == len(even_weight) == sum(dp[2 * n])
        assert len(set(odd_weight)) == len(odd_weight) == sum(dp[2 * n + 1])
        for weight, collection in [(2 * n, even_weight), (2 * n + 1, odd_weight)]:
            for parts in collection:
                assert sum(parts) == weight and all(x > 0 for x in parts)
                assert all(a >= b for a, b in zip(parts, parts[1:]))
                assert reverse_sum(parts) == recurrence_sum(parts)
                paired = all(parts[i] == parts[i + 1] for i in range(0, len(parts) - 1, 2))
                assert (reverse_sum(parts) >= 0) == (len(parts) % 2 == 1 or paired)
                if len(parts) % 2:
                    assert reverse_sum(parts) > 0
                if parts:
                    rhs = sum(x % 2 for x in conjugate(parts))
                    if len(parts) % 2 == 0:
                        rhs = -rhs
                    assert reverse_sum(parts) == rhs
        a = [parts for parts in even_weight if reverse_sum(parts) >= 0]
        b = {parts for parts in odd_weight if len(parts) % 2}
        c = [parts for parts in even_weight if reverse_sum(parts) < 0]
        d = [parts for parts in odd_weight
             if sum(x % 2 for x in parts) % 2 == 1
             and sum(x % 2 == 0 for x in parts) % 2 == 1]
        assert set(d) == {parts for parts in odd_weight if len(parts) % 2 == 0}
        assert len(b) == dp[2 * n + 1][1]
        assert len(d) == dp[2 * n + 1][0]
        images = Counter(injection(parts) for parts in a)
        assert all(image in b for image in images)
        assert all(multiplicity == 1 for multiplicity in images.values())
        odd_images = {injection(parts) for parts in a if len(parts) % 2}
        even_images = {injection(parts) for parts in a if len(parts) % 2 == 0}
        assert odd_images.isdisjoint(even_images)
        for parts in a:
            image = injection(parts)
            if len(parts) % 2:
                assert (image[0] - 1,) + image[1:] == parts
                assert len(image) == 1 or image[0] > image[1]
            else:
                assert image[-1] == 1 and image[:-1] == parts
                if n > 0:
                    assert image[0] == image[1]
        assert len(a) <= len(b)
        assert len(a) + len(c) == sum(dp[2 * n])
        assert len(b) + len(d) == sum(dp[2 * n + 1])
        if n >= 4:
            witness = (n - 1, n - 1, 3)
            assert witness in b and witness not in images
            assert len(a) < len(b)
        else:
            assert len(a) == len(b)
        if n < len(SOURCE_A):
            assert (len(a), len(b), len(c), len(d)) == (
                SOURCE_A[n], SOURCE_B[n], SOURCE_C[n], SOURCE_D[n])
        if n < len(SOURCE_DIFF):
            assert len(b) - len(a) == SOURCE_DIFF[n]
        if n == 4:
            assert set(a) == SOURCE_EXAMPLE_WEIGHT_8
        checked_sources += len(a)
        rows.append((n, len(a), len(b), len(c), len(d), sum(dp[2 * n]), sum(dp[2 * n + 1])))
    print('Exact enumeration: n=0..%d; source fidelity: four OEIS prefixes, A344611 difference and weight-8 example.' % MAX_N)
    print('n A344611 A160786 A344743 A236914 p(2n) p(2n+1)')
    for row in rows:
        print(*row)
    print('Map validity, inverse recovery, injectivity and branch disjointness:', checked_sources, 'source partitions checked.')
    print('Strictness witness (n-1,n-1,3): checked for n=4..%d; equality at n=0..3.' % MAX_N)
    assert rows[1][3] == 0 and rows[1][4] == 1
    assert rows[1][5] == 2 and rows[1][6] == 3
    print('Complement identity: A344611+A344743=p(2n), but A160786+A236914=p(2n+1).')
    print('A344743 conjecture is not equivalent: at n=1 it asserts 0 >= 1, which is false.')
    print('Finite checks only; no Lean proof of claim or its negation is supplied.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
