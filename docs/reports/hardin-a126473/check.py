#!/usr/bin/env python3
"""Exact finite checks for Hardin's two empirical A126473 alphabet increments.

Only Python integer arithmetic is used. These are finite checks, not a proof.
"""
from collections import Counter
from itertools import product

# OEIS A126473 revision 57, %S/%T/%U, indexed by word length starting at zero.
SOURCE_TERMS = (
    1, 5, 23, 107, 497, 2309, 10727, 49835, 231521, 1075589,
    4996919, 23214443, 107848529, 501037445, 2327695367, 10813893803,
    50238661313, 233396326661, 1084301290583, 5037394142315,
    23402480441009, 108722104190981, 505095858086951, 2346549744920747,
)


def count_dp(base, length, radius=3):
    """Count source words by their last letter, with leading zeros allowed."""
    assert base >= 0 and length >= 0 and radius >= 1
    if length == 0:
        return 1
    ending = [1] * base
    for _ in range(1, length):
        ending = [sum(ending[max(0, j - radius):min(base, j + radius + 1)])
                  for j in range(base)]
    return sum(ending)


def count_literal(base, length):
    """Independent enumeration of the source's original non-circular words."""
    return sum(all(abs(v - u) <= 3 for u, v in zip(word, word[1:]))
               for word in product(range(base), repeat=length))


def spans(length, radius):
    """Enumerate increment words; include the empty prefix at position zero."""
    histogram = Counter()
    maximal_words = []
    for increments in product(range(-radius, radius + 1), repeat=length - 1):
        position = low = high = 0
        for step in increments:
            position += step
            low = min(low, position)
            high = max(high, position)
        span = high - low
        assert 0 <= span <= radius * (length - 1)
        histogram[span] += 1
        if span == radius * (length - 1):
            maximal_words.append(increments)
    if length >= 2:
        assert set(maximal_words) == {
            (-radius,) * (length - 1), (radius,) * (length - 1)}
        assert len(maximal_words) == 2
    else:
        assert maximal_words == [()]
    assert sum(histogram.values()) == (2 * radius + 1) ** (length - 1)
    return histogram


def main():
    actual = tuple(count_dp(5, n) for n in range(len(SOURCE_TERMS)))
    assert actual == SOURCE_TERMS
    assert all(actual[n] == 4 * actual[n - 1] + 3 * actual[n - 2]
               for n in range(2, len(actual)))
    print('FIDELITY: A(5,n) matches all 24 source terms, n=0..23')
    print('FIDELITY_PREFIX:', ','.join(map(str, actual[:9])))

    literal_cases = 0
    for n in range(1, 6):
        for b in range(0, 9):
            assert count_literal(b, n) == count_dp(b, n)
            literal_cases += 1
    print(f'LITERAL_WORDS: {literal_cases} cases, n=1..5 and b=0..8, agree with DP')

    first = second = 0
    for n in range(1, 13):
        for b in range(1, 3 * n + 5):
            value = count_dp(b, n)
            previous = count_dp(b - 1, n)
            if b >= 3 * n - 2:
                assert value == previous + 7 ** (n - 1)
                first += 1
            if n >= 2 and b == 3 * n - 3:
                assert value == previous + 7 ** (n - 1) - 2
                second += 1
    print(f'PUBLISHED_CLAUSES: first={first}, second={second}; n=1..12, b=1..3n+4')

    fiber_cases = 0
    for n in range(1, 8):
        histogram = spans(n, 3)
        for b in range(0, 3 * n + 5):
            assert count_dp(b, n) == sum(
                multiplicity * max(b - span, 0)
                for span, multiplicity in histogram.items())
            if b >= 1:
                assert count_dp(b, n) - count_dp(b - 1, n) == sum(
                    multiplicity for span, multiplicity in histogram.items()
                    if span <= b - 1)
            fiber_cases += 1
    print(f'INCREMENT_ENUMERATION: {fiber_cases} fiber cases, n=1..7; extremal spans checked')

    general_cases = 0
    for radius in range(1, 5):
        for n in range(1, 11):
            threshold = radius * (n - 1)
            for b in range(1, threshold + 7):
                delta = count_dp(b, n, radius) - count_dp(b - 1, n, radius)
                if b >= threshold + 1:
                    assert delta == (2 * radius + 1) ** (n - 1)
                if n >= 2 and b == threshold:
                    assert delta == (2 * radius + 1) ** (n - 1) - 2
                general_cases += 1
    print(f'GENERAL_STEP_BOUND: {general_cases} DP cases, radius=1..4, n=1..10')
    for n, b, expected in [(2, 3, 5), (2, 4, 7), (3, 6, 47), (3, 7, 49)]:
        high, low = count_dp(b, n), count_dp(b - 1, n)
        assert high - low == expected
        print(f'BOUNDARY: n={n}, b={b}, A(b,n)={high}, A(b-1,n)={low}, delta={expected}')
    assert all(count_dp(0, n) == 0 for n in range(1, 13))
    print('DEGENERATE_ALPHABET: A(0,n)=0 for n=1..12; length zero excluded from claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
