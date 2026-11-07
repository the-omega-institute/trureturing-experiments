#!/usr/bin/env python3
"""Exact bounded checks for OEIS A183625's native arrays and Barker formulas.

Source: https://oeis.org/A183625/internal, revision 10.
No recurrence or guessed formula is used to generate the native counts.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


ROWS = tuple(product(range(3), repeat=3))


def compatible(r, s):
    return all(r[j] + r[j + 1] + s[j] + s[j + 1] == 4 for j in range(2))


def admissible(array):
    return all(compatible(r, s) for r, s in zip(array, array[1:]))


def profile(row):
    return (row[0] + row[1], row[1] + row[2])


def complement_odd(array):
    return tuple(tuple(2 - x for x in r) if i % 2 else r
                 for i, r in enumerate(array))


def closed_form(n):
    assert n >= 1
    return 3 * (4 + 2 ** (2 + n) + 3 ** n)


def transfer_counts(max_height):
    # The state is the entire ordered row. Test the original four entries.
    next_rows = {r: tuple(s for s in ROWS if compatible(r, s)) for r in ROWS}
    weights = {r: 1 for r in ROWS}
    counts = [0, len(ROWS)]
    for _height in range(2, max_height + 1):
        new = {r: 0 for r in ROWS}
        for r, weight in weights.items():
            for s in next_rows[r]:
                new[s] += weight
        weights = new
        counts.append(sum(weights.values()))
    return counts


def poly_mul(p, q):
    out = [Fraction(0) for _ in range(len(p) + len(q) - 1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def main():
    source = {}
    for line in Path(__file__).with_name('source-terms.txt').read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            n, value = map(int, line.split())
            source[n] = value
    assert sorted(source) == list(range(1, 201))
    printed_prefix = [45, 87, 189, 447, 1125, 2967, 8109, 22767, 65205,
                      189447, 556029, 1643487]
    counts = transfer_counts(201)
    assert counts[2:14] == printed_prefix
    for n, value in source.items():
        assert counts[n + 1] == value == closed_form(n), (n, value, counts[n + 1])
    print('SOURCE_FIDELITY: native transfer counts reproduce the 12 printed initial terms')
    print('SOURCE_BFILE: all 200 source terms n=1..200 match native transfer counts and closed form')
    print('HEIGHTS_1_TO_12:', ','.join(map(str, counts[1:13])))

    exhaustive = []
    for h in range(1, 5):
        total = sum(1 for array in product(ROWS, repeat=h) if admissible(array))
        assert total == counts[h] == 12 + 6 * 2 ** h + 3 ** h
        exhaustive.append(total)
    print('EXHAUSTIVE_NATIVE: heights 1..4, ambient arrays up to 531441; counts=' + str(exhaustive))

    fibers = defaultdict(list)
    for r in ROWS:
        fibers[profile(r)].append(r)
    table = [[len(fibers[(s, t)]) for t in range(5)] for s in range(5)]
    expected_table = [[1, 1, 1, 0, 0], [1, 2, 2, 1, 0], [1, 2, 3, 2, 1],
                      [0, 1, 2, 2, 1], [0, 0, 1, 1, 1]]
    assert table == expected_table
    for s in range(5):
        for t in range(5):
            assert table[s][t] == max(0, min(2, s, t) - max(0, s - 2, t - 2) + 1)
    assert Counter(map(len, fibers.values())) == Counter({1: 12, 2: 6, 3: 1, 0: 6})
    for r, s in product(ROWS, repeat=2):
        pr, ps = profile(r), profile(s)
        assert compatible(r, s) == (pr[0] + ps[0] == 4 and pr[1] + ps[1] == 4)
    print('ROW_FIBERS:', table)
    print('PROFILE_TRANSITION: all 729 ordered row pairs checked against native block equations')

    for h in range(1, 4):
        images = set()
        for array in product(ROWS, repeat=h):
            if admissible(array):
                transformed = complement_odd(array)
                p = profile(array[0])
                assert all(profile(r) == p for r in transformed)
                assert complement_odd(transformed) == array
                images.add((p, transformed))
        targets = {(p, rows) for p, fiber in fibers.items()
                   for rows in product(fiber, repeat=h)}
        assert images == targets
        for _p, rows in targets:
            assert admissible(complement_odd(rows))
        assert len(images) == counts[h]
    print('BIJECTION_CHECK: both directions, heights 1..3, with labelled positions preserved')

    # Formal-series equality is checked coefficientwise with exact Fractions.
    denominator = poly_mul(poly_mul([1, -1], [1, -2]), [1, -3])
    numerator = [0, 45, -183, 162]
    assert denominator == [1, -6, 11, -6]
    coefficients = [Fraction(0)] + [Fraction(counts[n + 1]) for n in range(1, 201)]
    for n in range(201):
        coefficient = sum(denominator[j] * coefficients[n - j]
                          for j in range(min(n, 3) + 1))
        assert coefficient == (numerator[n] if n < len(numerator) else 0)
    partial_numerator = [Fraction(0) for _ in range(4)]
    for weight, others in [(12, (2, 3)), (24, (1, 3)), (9, (1, 2))]:
        term = poly_mul([0, weight], poly_mul([1, -others[0]], [1, -others[1]]))
        partial_numerator = [a + b for a, b in zip(partial_numerator, term)]
    assert partial_numerator == numerator
    assert denominator[0] == 1
    print('GENERATING_FUNCTION: exact product coefficients degrees 0..200; numerator [0,45,-183,162]')
    for n in range(4, 201):
        assert counts[n + 1] == 6 * counts[n] - 11 * counts[n - 1] + 6 * counts[n - 2]
    print('EMPIRICAL_RECURRENCE: reproduced for n=4..200 from native transfer counts')

    for w in range(1, 9):
        groups = Counter(tuple(row[j] + row[j + 1] for j in range(w - 1))
                         for row in product(range(3), repeat=w))
        multiplicities = Counter(groups.values())
        expected = Counter({1: 3 ** w - 2 ** (w + 1) + 1,
                            2: 2 ** w - 2, 3: 1})
        assert +multiplicities == +expected, (w, multiplicities, expected)
    print('EXTENSION_CHECK: profile-fiber multiplicities for widths 1..8, exact enumeration')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
