#!/usr/bin/env python3
"""Exact source-fidelity and ratio checks; this is not an all-parameter proof."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb


def forms(n, k):
    """Section 5, equation (2), in the label order I, II, III, IV, V, VI."""
    return tuple((a % n, b % n) for a, b in
                 ((1, k - 1), (k, 1 - k), (-1, k),
                  (k - 1, -k), (-k, 1), (1 - k, -1)))


def allowed(n, k):
    pairs = forms(n, k)
    return tuple(pairs[i] for i in ((1, 2, 3, 4) if k == 2 else range(6)))


def histogram(n, classes):
    """Count labelled partial choices by the full integer imbalance vector."""
    counts = Counter()
    imbalance = [0] * n

    def visit(i):
        if i == len(classes):
            counts[tuple(imbalance[1:])] += 1
            return
        for a, b in allowed(n, classes[i]):
            assert a != 0 and b != 0
            imbalance[a] += 1
            imbalance[b] -= 1
            visit(i + 1)
            imbalance[a] -= 1
            imbalance[b] += 1

    visit(0)
    return counts


def balanced_mitm(t):
    """Join opposite imbalance vectors, with no closure/connectivity test."""
    assert t >= 2
    n = 3 * t
    classes = tuple(range(2, n - 1))
    split = len(classes) // 2
    left = histogram(n, classes[:split])
    right = histogram(n, classes[split:])
    return sum(multiplicity * right[tuple(-v for v in vector)]
               for vector, multiplicity in left.items())


def balanced_direct(t):
    """Separate whole-choice enumeration for the two smallest moduli."""
    n = 3 * t
    answer = 0
    for plan in product(*(allowed(n, k) for k in range(2, n - 1))):
        outgoing = Counter(a for a, _ in plan)
        incoming = Counter(b for _, b in plan)
        answer += outgoing == incoming
    return answer


def closed_count(t):
    """The frozen repository closed forms, not a new definition of the count."""
    assert t >= 2
    r = t // 2
    if t % 2 == 0:
        return 4 * 6 ** (r - 1) * comb(2 * r - 2, r - 1)
    return 2 * 6 ** r * comb(2 * r - 1, r)


def main():
    source_table = (4, 12, 48, 216, 864, 4320, 17280, 90720,
                    362880, 1959552, 7838208)
    computed = tuple(closed_count(t) for t in range(2, 13))
    assert computed == source_table
    print('Source Conjecture 9.3 table, t=2..12:', ', '.join(map(str, computed)))
    for t in range(2, 6):
        count = balanced_mitm(t)
        assert count == source_table[t - 2] == closed_count(t)
        if t <= 3:
            assert balanced_direct(t) == count
        print(f'Independent labelled balance count: t={t}, n={3*t}, count={count}')
    for r in range(1, 201):
        even = Fraction(closed_count(2 * r + 2), closed_count(2 * r))
        odd = Fraction(closed_count(2 * r + 3), closed_count(2 * r + 1))
        assert even == Fraction(24 * (2 * r - 1), 2 * r)
        assert odd == Fraction(12 * (2 * r + 1), r + 1)
        assert 24 - even == Fraction(12, r)
        assert 24 - odd == Fraction(12, r + 1)
        assert 0 < even < 24 and 0 < odd < 24
        assert 0 < closed_count(2 * r) and 0 < closed_count(2 * r + 1)
    print('Exact even and odd ratio/error identities: r=1..200, all passed')
    for t in range(2, 401):
        ratio = Fraction(closed_count(t + 2), closed_count(t))
        assert abs(ratio - 24) <= Fraction(24, t)
    print('Uniform exact bound |a(t+2)/a(t)-24| <= 24/t: t=2..400, all passed')
    assert Fraction(24, 27) == Fraction(8, 9)
    print('Source charged-sector normalization: 24/27 = 8/9')
    print('Finite checks support the outline; the unbounded limit is not certified here.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
