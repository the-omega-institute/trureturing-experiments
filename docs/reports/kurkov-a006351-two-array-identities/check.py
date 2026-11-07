#!/usr/bin/env python3
"""Exact fidelity and finite identity checks for OEIS A006351/A379459/A379460.

All calculations use integers or fractions.Fraction. This is a finite check,
not an unbounded proof. The recurrences are transcribed from the OEIS NAME
and FORMULA fields; no source program is copied.
"""
from fractions import Fraction
from math import comb, factorial

LIMIT = 30
DEGREE = LIMIT + 2
WIDTH = LIMIT + 9

# Published DATA, in downward antidiagonal order (first 28 entries).
SOURCE_U = [2, 4, 8, 12, 36, 52, 48, 192, 368, 472, 240, 1200,
            2880, 4560, 5504, 1440, 8640, 24960, 47280, 67408,
            78416, 10080, 70560, 238560, 527520, 871584, 1163232, 1320064]
SOURCE_V = [1, 1, 2, 1, 6, 8, 1, 12, 44, 52, 1, 20, 140, 420,
            472, 1, 30, 340, 1860, 5032, 5504, 1, 42, 700, 6020,
            28672, 72912, 78416]
SOURCE_A = [1, 2, 8, 52, 472, 5504, 78416, 1320064, 25637824,
            564275648, 13879795712]


def arrays():
    u = [[Fraction(2 * factorial(k + 1)) for k in range(WIDTH)]]
    v = [[1] * WIDTH]
    for _ in range(LIMIT):
        prev = u[-1]
        u.append([prev[k + 1] + 2 * factorial(k + 1) *
                  sum((prev[j] / factorial(j) for j in range(k + 1)), Fraction(0))
                  for k in range(len(prev) - 1)])
        prev_v = v[-1]
        row = [prev_v[0] + prev_v[1]]
        for k in range(1, len(prev_v) - 1):
            row.append(row[k - 1] + (k + 1) * (prev_v[k] + prev_v[k + 1]))
        v.append(row)
    assert all(x.denominator == 1 for row in u for x in row)
    return u, v


def sequence():
    # OEIS A006351, Ilya Gutkovskiy's independent FORMULA recurrence.
    a = [0, 1]
    for n in range(2, DEGREE + 1):
        a.append(a[n - 1] + sum(comb(n - 1, k) * a[k] * a[n - k]
                                for k in range(1, n)))
    return a


def product(p, q):
    return [sum((p[j] * q[n - j] for j in range(n + 1)), Fraction(0))
            for n in range(DEGREE + 1)]


def check_reversion(a):
    f = [Fraction(a[n], factorial(n)) for n in range(DEGREE + 1)]
    # Compose 2 log(1+x)-x with F, coefficientwise modulo t^(DEGREE+1).
    power = [Fraction(1)] + [Fraction(0)] * DEGREE
    result = [Fraction(0)] * (DEGREE + 1)
    for j in range(1, DEGREE + 1):
        power = product(power, f)
        c = 2 * Fraction((-1) ** (j + 1), j) - (1 if j == 1 else 0)
        result = [old + c * x for old, x in zip(result, power)]
    assert result == [Fraction(0), Fraction(1)] + [Fraction(0)] * (DEGREE - 1)


def antidiagonals(rows):
    return [rows[r][d - r] for d in range(7) for r in range(d + 1)]


def main():
    u, v = arrays()
    a = sequence()
    assert antidiagonals(u) == SOURCE_U
    assert antidiagonals(v) == SOURCE_V
    assert a[1:len(SOURCE_A) + 1] == SOURCE_A
    print('FIDELITY: 28 A379459 terms, 28 A379460 terms, 11 A006351 terms match OEIS DATA')
    check_reversion(a)
    print(f'EGF: (2 log(1+x)-x)(F(t)) = t through degree {DEGREE}, exactly over Q')
    for r in range(LIMIT + 1):
        assert u[r][0] == a[r + 2], ('U', r)
        assert v[r][0] == a[r + 1], ('V', r)
    print(f'IDENTITIES: U(r,0)=a(r+2) and V(r,0)=a(r+1) for every r=0..{LIMIT}')
    print('U(r,0), r=0..6:', ', '.join(str(int(u[r][0])) for r in range(7)))
    print('V(r,0), r=0..6:', ', '.join(str(v[r][0]) for r in range(7)))
    print('a(32) =', a[32])
    print('SCOPE: finite indexing/fidelity checks; no proof of the unbounded conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
