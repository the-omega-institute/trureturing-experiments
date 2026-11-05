#!/usr/bin/env python3
"""Exact independent checks for Barker's OEIS A097399 generating function.

Only Python's standard library is used. Integers and Fractions carry all evidence.
The finite checks support a proposed proof; they are not a Lean proof.
"""

from fractions import Fraction
from itertools import permutations
from math import factorial


VALUES = tuple(range(-4, 5))
ENDPOINTS = (0, 1, 2, 12)
BOUNDS = (86, 104, 172, 972, 81)
SOURCE_TERMS = (
    86, 104, 172, 252, 332, 412, 492, 572, 652, 732, 812, 892, 972,
    1053, 1134, 1215, 1296, 1377, 1458, 1539, 1620, 1701, 1782,
    1863, 1944, 2025, 2106, 2187, 2268, 2349, 2430, 2511, 2592,
    2673, 2754, 2835, 2916, 2997, 3078, 3159, 3240, 3321, 3402,
    3483, 3564,
)
WITNESSES = (
    (-4, -3, 0, 1, -1, 4, -2, 3, 2),
    (-4, -3, 1, 2, -1, 3, -2, 4, 0),
    (-4, -1, 3, 2, -3, 1, 0, 4, -2),
    (-4, -1, 2, 3, -3, 0, 1, 4, -2),
)
WITNESS_POLYS = ((86, 4, 0, 0), (72, 32, 0, 0),
                 (12, 80, 0, 0), (0, 81, 0, 0))


def determinant(b):
    a, c, e, f, g, h, i, j, k = b
    return a * (g * k - h * j) - c * (f * k - h * i) + e * (f * j - g * i)


def shifted(b, n):
    return tuple(x + n for x in b)


def determinant_gauss(b):
    """Independent exact rational Gaussian elimination, with row pivoting."""
    a = [[Fraction(x) for x in b[i:i + 3]] for i in (0, 3, 6)]
    result = Fraction(1)
    for col in range(3):
        pivot = next((r for r in range(col, 3) if a[r][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            result = -result
        q = a[col][col]
        result *= q
        for row in range(col + 1, 3):
            ratio = a[row][col] / q
            for j in range(col + 1, 3):
                a[row][j] -= ratio * a[col][j]
            a[row][col] = Fraction(0)
    return result


def poly_product(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def determinant_poly(b):
    """Compute every coefficient of det(B+tJ) by the Leibniz definition."""
    result = [0, 0, 0, 0]
    for p in permutations(range(3)):
        inversions = sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
        term = [1]
        for i in range(3):
            term = poly_product(term, [b[3 * i + p[i]], 1])
        for i, x in enumerate(term):
            result[i] += (-1 if inversions % 2 else 1) * x
    return tuple(result)


def normalized(b):
    return b[0] == -4 and b[1] < b[2] and b[3] < b[6] and b[1] < b[3]


def normalize(b):
    m = [list(b[i:i + 3]) for i in (0, 3, 6)]
    r, c = divmod(b.index(-4), 3)
    m[0], m[r] = m[r], m[0]
    for row in m:
        row[0], row[c] = row[c], row[0]
    if m[0][1] > m[0][2]:
        for row in m:
            row[1], row[2] = row[2], row[1]
    if m[1][0] > m[2][0]:
        m[1], m[2] = m[2], m[1]
    if m[0][1] > m[1][0]:
        m = list(zip(*m))
    out = tuple(x for row in m for x in row)
    assert normalized(out)
    return out


def predicted(n):
    if n == 0:
        return 86
    if n == 1:
        return 104
    return 80 * n + 12 if n <= 12 else 81 * n


def main():
    # Direct unrestricted enumeration independently checks the five finite bounds.
    maxima = [0] * 5
    signed_maxima = [None] * 4
    signed_minima = [None] * 4
    representatives = set()
    count = 0
    for b in permutations(VALUES):
        count += 1
        vals = tuple(determinant(shifted(b, n)) for n in ENDPOINTS)
        slope = vals[1] - vals[0]
        values = tuple(abs(x) for x in vals) + (abs(slope),)
        maxima = [max(old, new) for old, new in zip(maxima, values)]
        signed_maxima = [v if old is None else max(old, v)
                         for old, v in zip(signed_maxima, vals)]
        signed_minima = [v if old is None else min(old, v)
                         for old, v in zip(signed_minima, vals)]
        rep = normalize(b)
        representatives.add(rep)
        assert tuple(abs(determinant(shifted(rep, n))) for n in ENDPOINTS) == values[:4]
        assert abs(determinant(shifted(rep, 1)) - determinant(rep)) == abs(slope)
    assert count == factorial(9)
    assert tuple(maxima) == BOUNDS
    assert signed_maxima == list(BOUNDS[:4])
    assert signed_minima == [-x for x in BOUNDS[:4]]
    print('Unrestricted arrangements:', count)
    print('Endpoint/slope absolute maxima (t=0,1,2,12; slope):', maxima)
    print('Signed endpoint maxima:', signed_maxima)
    print('Signed endpoint minima:', signed_minima)

    reduced = {(-4,) + b for b in permutations(VALUES[1:]) if normalized((-4,) + b)}
    assert len(reduced) == 5040 and representatives == reduced
    gauss_maxima = [0] * 5
    small_maxima = [0] * 46
    for b in sorted(reduced):
        c, s, quadratic, cubic = determinant_poly(b)
        assert quadratic == cubic == 0
        assert c == determinant(b)
        assert s == determinant(shifted(b, 1)) - c
        vals = tuple(determinant_gauss(shifted(b, n)) for n in ENDPOINTS)
        assert vals == tuple(c + n * s for n in ENDPOINTS)
        values = tuple(abs(x) for x in vals) + (abs(vals[1] - vals[0]),)
        gauss_maxima = [max(old, new) for old, new in zip(gauss_maxima, values)]
        for n in range(46):
            actual = abs(determinant(shifted(b, n)))
            assert actual == abs(c + s * n)
            small_maxima[n] = max(small_maxima[n], actual)
    assert tuple(gauss_maxima) == BOUNDS
    assert tuple(small_maxima[:45]) == SOURCE_TERMS
    assert small_maxima == [predicted(n) for n in range(46)]
    print('Normalized representatives:', len(reduced))
    print('Rational Gaussian endpoint/slope maxima:', [int(x) for x in gauss_maxima])
    print('All normalized determinant polynomials have quadratic=cubic=0')
    print('Fidelity: all 45 OEIS revision-21 terms reproduced, including a(0)=86 and a(5)=412')
    print('Direct small-shift maxima n=0..45:', small_maxima)

    for b, p in zip(WITNESSES, WITNESS_POLYS):
        assert tuple(sorted(b)) == VALUES
        assert determinant_poly(b) == p
    print('Attaining witness polynomials (constant,slope,quadratic,cubic):', WITNESS_POLYS)
    assert determinant_gauss(WITNESSES[0]) == 86
    assert determinant_gauss(shifted(WITNESSES[2], 5)) == 412

    # Exact coefficient identity on a larger window, including both transitions.
    second_difference = {}
    for n in range(1001):
        value = predicted(n)
        if n >= 1:
            value -= 2 * predicted(n - 1)
        if n >= 2:
            value += predicted(n - 2)
        if value:
            second_difference[n] = value
    assert second_difference == {0: 86, 1: -68, 2: 50, 3: 12, 13: 1}
    print('Nonzero coefficients of (1-x)^2*A through degree 1000:', second_difference)
    print('Scope: exact finite checks; the unbounded proof is a preregistered outline, not Lean-verified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
