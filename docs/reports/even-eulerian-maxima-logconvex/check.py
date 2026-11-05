#!/usr/bin/env python3
"""Exact checks for Ambrus--Gárgyán's even normalized Eulerian-maxima conjecture."""
from fractions import Fraction
from itertools import permutations
from math import factorial


def ascents(p):
    return sum(a < b for a, b in zip(p, p[1:]))


def row_by_enumeration(m):
    row = [0] * m
    for p in permutations(range(1, m + 1)):
        row[ascents(p)] += 1
    return row


def eulerian_row(m):
    # A(m,l) counts m-permutations with l-1 adjacent strict ascents.
    row = [1]
    for n in range(1, m + 1):
        nxt = [0] * n
        for asc, count in enumerate(row):
            # Insert the new maximum into n positions: one extra ascent
            # in n-1-asc positions and no extra ascent in asc+1 positions.
            if asc < n:
                nxt[asc] += (asc + 1) * count
            if asc + 1 < n:
                nxt[asc + 1] += (n - asc - 1) * count
        row = nxt
    return row


def f(m):
    row = eulerian_row(m)
    return Fraction(max(row), factorial(m))

# Fidelity: direct permutation counts agree with the recurrence for small rows.
for m in range(1, 9):
    assert row_by_enumeration(m) == eulerian_row(m), (m, row_by_enumeration(m), eulerian_row(m))

# The source's standard Eulerian maxima and the first conjectural center.
assert max(eulerian_row(4)) == 11
assert max(eulerian_row(6)) == 302
assert max(eulerian_row(8)) == 15619
assert f(4) == Fraction(11, 24)
assert f(6) == Fraction(151, 360)
assert f(8) == Fraction(15619, 40320)
assert f(4) * f(8) - f(6) ** 2 == Fraction(23423, 14515200)
assert f(2) * f(6) - f(4) ** 2 == Fraction(-1, 2880)

# Check every centered comparison through p=100 with exact rational arithmetic.
for p in range(3, 101):
    assert f(2 * p) ** 2 <= f(2 * p - 2) * f(2 * p + 2), p

# Check the published ratio envelope and the strict separation used in the outline.
for p in range(4, 101):
    ratio = f(2 * p) / f(2 * p - 2)
    lower = Fraction(2 * p * p, 2 * p * p + p - 1)
    upper = Fraction(p * (2 * p + 1) ** 2, 2 * (p + 1) ** 2 * (2 * p - 1))
    assert lower < ratio < upper, (p, ratio, lower, upper)
    numerator = 2 * p ** 3 - 5 * p ** 2 - 10 * p - 4
    t = p - 4
    assert numerator == 2 * t ** 3 + 19 * t ** 2 + 46 * t + 4
    assert numerator > 0

# The source's weaker product lower bound is checked at the same indices.
for p in range(3, 101):
    product = f(2 * p + 2) / f(2 * p) * f(2 * p - 2) / f(2 * p)
    bound = Fraction(4 * p ** 3 + 6 * p ** 2 - 2, 4 * p ** 3 + 6 * p ** 2 - 1)
    assert product >= bound, p

print("checked direct permutations m=1..8 and recurrence through m=202")
print("checked f(4), f(6), f(8), the sharp p=2 failure, ratio envelopes, and polynomial identity")
print("ALL_CHECKS_OK")
