#!/usr/bin/env python3
"""Exact preregistration checks for Shao--Zhao arXiv:2610.00108v1, section 9.

All arithmetic used as evidence is integral or Fraction arithmetic. No packages
outside the Python 3 standard library are required. The proposed formula counts
labelled functions, with no quotient by rotations or reflections.
"""

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path


def C(top, bottom):
    """Ordinary binomial coefficient, zero outside 0 <= bottom <= top."""
    return comb(top, bottom) if 0 <= bottom <= top else 0


def closed_form(n, j):
    """Proposed finite sum, evaluated without recursion or floating point."""
    assert n >= 3 and 0 <= j <= 2 * n
    r, s = j - n, 2 * n - j
    zero_free = Fraction(n * C(s, r), s) if 0 <= 2 * r <= n else Fraction(0)
    exceptional = 3 if n % 3 == 0 and 3 * j == 4 * n else 0
    answer = zero_free + exceptional
    for b in range(j // 4 + 1):
        for a in range((j - 4 * b) // 2 + 1):
            c = j - 2 * a - 4 * b
            t = n - j + a + b
            if t < 1:
                continue
            m = t + b
            assert m > 0 and c >= 0
            term = Fraction(n * C(m, b) * C(c + t - 1, c) * C(c + t, a), m)
            assert term.denominator == 1, (n, j, a, b, term)
            answer += term
    assert answer.denominator == 1, (n, j, answer)
    return answer.numerator


def polynomial_form(n, j):
    """Optional fixed-deficit corollary in its stable range n > j."""
    assert n >= 3 and 0 <= j < n
    answer = Fraction(C(n, j))
    for a in range(1, j // 2 + 1):
        answer += Fraction(n * C(n - a, j - 2 * a) * C(n - a - 1, a - 1), a)
    for b in range(1, j // 4 + 1):
        for a in range((j - 4 * b) // 2 + 1):
            c = j - 2 * a - 4 * b
            answer += Fraction(
                n * C(n - j + a + 2 * b - 1, b - 1)
                * C(n - a - 3 * b - 1, c) * C(n - a - 3 * b, a), b)
    assert answer.denominator == 1
    return answer.numerator


def italian(f):
    n = len(f)
    return all(f[i] != 0 or f[(i - 1) % n] + f[(i + 1) % n] >= 2
               for i in range(n))


def brute_row(n):
    """Directly use f and the source's neighbourhood inequality."""
    row = [0] * (2 * n + 1)
    for f in product(range(3), repeat=n):
        if italian(f):
            row[sum(f)] += 1
    return row


def weighted_dp_row(n):
    """Independent boundary-conditioned DP in the original f variables.

    State (a,b,k) retains the two latest labels and the weight. Adding c checks
    the middle vertex b; the two remaining cyclic conditions are checked last.
    """
    row = [0] * (2 * n + 1)
    for first, second in product(range(3), repeat=2):
        states = {(first, second, first + second): 1}
        for _ in range(2, n):
            nxt = defaultdict(int)
            for (a, b, k), count in states.items():
                for c in range(3):
                    if b != 0 or a + c >= 2:
                        nxt[b, c, k + c] += count
            states = nxt
        for (a, b, k), count in states.items():
            if (b != 0 or a + first >= 2) and (first != 0 or b + second >= 2):
                row[k] += count
    return row


# Source section 8.3, Tables 4--5 (labels tab:cycle_a/tab:cycle_b).
# Index k is the original weight. Omitted coefficients, including k=0, are zero.
SOURCE_TABLES = {
    3: [0, 0, 6, 7, 6, 3, 1],
    4: [0, 0, 2, 8, 19, 16, 10, 4, 1],
    5: [0, 0, 0, 5, 25, 41, 45, 30, 15, 5, 1],
    6: [0, 0, 0, 2, 18, 60, 105, 114, 90, 50, 21, 6, 1],
    7: [0, 0, 0, 0, 7, 49, 154, 260, 301, 252, 161, 77, 28, 7, 1],
    8: [0, 0, 0, 0, 2, 24, 144, 384, 659, 784, 704, 488, 266, 112, 36, 8, 1],
    9: [0, 0, 0, 0, 0, 9, 87, 387, 990, 1663, 2052, 1935, 1446, 864, 414, 156, 45, 9, 1],
    10: [0, 0, 0, 0, 0, 2, 35, 270, 1055, 2530, 4253, 5350, 5285, 4190, 2710, 1432, 615, 210, 55, 10, 1],
    11: [0, 0, 0, 0, 0, 0, 11, 132, 814, 2827, 6523, 10891, 13992, 14322, 11968, 8261, 4741, 2255, 880, 275, 66, 11, 1],
    12: [0, 0, 0, 0, 0, 0, 2, 48, 462, 2360, 7572, 16836, 28033, 36588, 38670, 33760, 24666, 15168, 7858, 3408, 1221, 352, 78, 12, 1],
}


def source_high_weight(n, j):
    """Source Theorem 7.1, j=0,...,5; the j=5 formula assumes n>=4."""
    unrestricted = sum(C(n, r) * C(n - r, j - 2 * r) for r in range(j // 2 + 1))
    return unrestricted - (2 * n if j == 5 else 0)


def main():
    rows = {n: weighted_dp_row(n) for n in range(3, 61)}
    for n, expected in SOURCE_TABLES.items():
        assert rows[n] == expected, ('source table', n, rows[n], expected)
    fidelity_count = 0
    for n in range(3, 61):
        for j in range(6):
            if j == 5 and n == 3:
                continue
            assert rows[n][2 * n - j] == source_high_weight(n, j), ('Theorem 7.1', n, j)
            fidelity_count += 1
    print('FIDELITY: source cycle tables n=3..12 reproduced in full')
    print(f'FIDELITY: Theorem 7.1 j=0..5, source domain respected, {fidelity_count} coefficients')
    print('SOURCE_TOTALS n=3..12:', [sum(rows[n]) for n in range(3, 13)])

    words = 0
    for n in range(3, 11):
        brute = brute_row(n)
        assert rows[n] == brute, ('DP vs brute', n)
        assert brute == [closed_form(n, 2 * n - k) for k in range(2 * n + 1)], ('formula vs brute', n)
        words += 3 ** n
    print(f'BRUTE_FORCE: all {words} labelled ternary functions for n=3..10, every coefficient agrees')
    print('DEFICIT_ROW n=3:', list(reversed(rows[3])))
    print('DEFICIT_ROW n=4:', list(reversed(rows[4])))

    tested, polynomial_tested = 0, 0
    totals = {}
    for n in range(3, 61):
        formula_row = [closed_form(n, j) for j in range(2 * n + 1)]
        assert formula_row == list(reversed(rows[n])), ('all-n formula vs DP', n)
        totals[n] = sum(formula_row)
        tested += len(formula_row)
        for j in range(n):
            assert polynomial_form(n, j) == formula_row[j], ('polynomial form', n, j)
            polynomial_tested += 1
    print(f'FORMULA: n=3..60, 0<=j<=2n, {tested} coefficients agree with independent weighted DP')
    print(f'POLYNOMIAL_COROLLARY: n=3..60, 0<=j<n, {polynomial_tested} exact checks')

    assert [totals[n] for n in range(3, 8)] == [23, 60, 167, 467, 1297]
    for n in range(3, 56):
        assert totals[n + 5] == (2 * totals[n + 4] + 2 * totals[n + 3]
                                + totals[n + 2] - totals[n + 1] - totals[n]), n
    print('REPOSITORY_RECURRENCE: five initial values and all 53 recurrences with base n=3..55 agree')
    print('LIMIT: these finite exact checks do not establish a universal theorem; no Lean proof was run')
    print('CHECK_SHA256:', sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
