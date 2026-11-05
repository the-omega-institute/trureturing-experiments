#!/usr/bin/env python3
"""Exact finite checks for five Kotesovec distinct-partition conjectures.

The finite tests do not prove any limit. Integer/rational interval arithmetic
certifies logarithms and square roots; no binary floating point is used.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, isqrt, prod

SOURCE_PREFIXES = {'A022629': [1,
             1,
             2,
             5,
             7,
             15,
             25,
             43,
             64,
             120,
             186,
             288,
             463,
             695,
             1105,
             1728,
             2525,
             3741,
             5775,
             8244,
             12447,
             18302,
             26424,
             37827,
             54729,
             78330,
             111184,
             159538,
             225624,
             315415,
             444708,
             618666,
             858165,
             1199701,
             1646076,
             2288961,
             3150951,
             4303995,
             5870539,
             8032571,
             10881794,
             14749051,
             19992626],
 'A092484': [1,
             1,
             4,
             13,
             25,
             77,
             161,
             393,
             726,
             2010,
             3850,
             7874,
             16791,
             31627,
             69695,
             139560,
             255997,
             482277,
             986021,
             1716430,
             3544299,
             6507128,
             11887340,
             21137849,
             38636535,
             70598032,
             123697772,
             233003286,
             412142276,
             711896765,
             1252360770],
 'A265840': [1,
             1,
             8,
             35,
             91,
             405,
             1069,
             3799,
             8686,
             36744,
             86310,
             235776,
             686329,
             1605779,
             5230579,
             13191702,
             30608501,
             73907925,
             206052723,
             433747560,
             1324608945,
             2995740974,
             6973434054,
             15364943439,
             35816669079,
             86662644756,
             184871083828,
             502089539734,
             1098571699830],
 'A265841': [1,
             1,
             16,
             97,
             337,
             2177,
             7313,
             38529,
             108594,
             717186,
             2053522,
             7527458,
             30757155,
             88042387,
             448973459,
             1390503396,
             4087546309,
             12699966117,
             49599776261,
             124699632310,
             608410782855,
             1651128186296,
             4862631132392,
             13170300313769,
             39285370060347,
             130999461143020],
 'A265842': [1,
             1,
             32,
             275,
             1267,
             11925,
             51445,
             406183,
             1406614,
             14690040,
             51144366,
             251885088,
             1481359033,
             5108404955,
             42614629915,
             158222158038,
             588574803125,
             2360755022421,
             13255325882835,
             39266011999104,
             325719196861377,
             1031732678138822,
             3791401325667894]}
ANUMS = tuple(SOURCE_PREFIXES)
NMAX = 180
ENUM_MAX = 45
SCALE = 10 ** 40
TERMS = 48


def coefficients(r, bound):
    """Multiply the source's factors, truncating after degree bound."""
    a = [1] + [0] * bound
    for k in range(1, bound + 1):
        for n in range(bound, k - 1, -1):
            a[n] += k ** r * a[n - k]
    return a


def distinct_parts(n, least=1):
    """Independent recursive enumeration, with strictly increasing parts."""
    if n == 0:
        yield ()
    for k in range(least, n + 1):
        for tail in distinct_parts(n - k, k + 1):
            yield (k,) + tail


def max_products_by_cardinality(bound):
    mmax = isqrt(2 * bound)
    table = [[0] * (mmax + 1) for _ in range(bound + 1)]
    table[0][0] = 1
    for k in range(1, bound + 1):
        for n in range(bound, k - 1, -1):
            for m in range(1, mmax + 1):
                table[n][m] = max(table[n][m], k * table[n - k][m - 1])
    return table


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[1], -x[0])


def mul(x, y):
    products = [a * b for a in x for b in y]
    return (min(products), max(products))


def const(x):
    return (Q(x), Q(x))


def times(c, x):
    return mul(const(c), x)


def log_unit(num, den):
    """Enclose log(num/den), for 1 <= num/den <= 2.

    log y = 2 sum_{j>=0} z^(2j+1)/(2j+1), z=(y-1)/(y+1).
    Each retained term is rounded downward on a fixed rational grid.
    The tail is <= 2/(3^(2T+1)*(2T+1)*(1-1/9)).
    """
    assert den <= num <= 2 * den
    if num == den:
        return const(0)
    zn, zd = num - den, num + den
    pn, pd, total = zn, zd, 0
    for j in range(TERMS):
        total += 2 * SCALE * pn // ((2 * j + 1) * pd)
        pn *= zn * zn
        pd *= zd * zd
    tail = Q(18, 8 * (2 * TERMS + 1) * 3 ** (2 * TERMS + 1))
    return (Q(total, SCALE), Q(total + TERMS, SCALE) + tail)


LOG2 = log_unit(2, 1)


@lru_cache(None)
def log_int(n):
    assert isinstance(n, int) and n > 0
    k = n.bit_length() - 1
    return add(times(k, LOG2), log_unit(n, 1 << k))


@lru_cache(None)
def sqrt_int(n):
    assert isinstance(n, int) and n >= 0
    lower = isqrt(n * SCALE * SCALE)
    upper = lower if lower * lower == n * SCALE * SCALE else lower + 1
    assert lower * lower <= n * SCALE * SCALE <= upper * upper
    return (Q(lower, SCALE), Q(upper, SCALE))


def triangular_witness(n):
    m = (isqrt(8 * n + 1) - 1) // 2
    t = m * (m + 1) // 2
    assert t <= n < (m + 1) * (m + 2) // 2
    d = n - t
    parts = tuple(range(1, m)) + (m + d,)
    assert len(set(parts)) == len(parts)
    assert min(parts) >= 1 and sum(parts) == n
    assert prod(parts) >= factorial(m)
    return m, parts


def main():
    arrays = [coefficients(r, NMAX) for r in range(6)]
    for r, anum in enumerate(ANUMS, 1):
        expected = SOURCE_PREFIXES[anum]
        assert arrays[r][:len(expected)] == expected, anum
        print(f'FIDELITY {anum} r={r}: all {len(expected)} displayed OEIS terms agree')
    assert arrays[1][6] == 6 + 5 + 8 + 6 == 25
    assert arrays[2][6] == 36 + 25 + 64 + 36 == 161
    print('FIDELITY: the two OEIS n=6 worked examples agree')

    total_partitions = 0
    for n in range(ENUM_MAX + 1):
        partitions = list(distinct_parts(n))
        total_partitions += len(partitions)
        assert len(partitions) == arrays[0][n]
        for r in range(1, 6):
            assert sum(prod(s) ** r for s in partitions) == arrays[r][n]
        for s in partitions:
            m = len(s)
            assert m * (m + 1) <= 2 * n
            if m:
                assert prod(s) * m ** m <= n ** m
    print(f'EXACT ENUMERATION: n=0..{ENUM_MAX}, {total_partitions} partitions; all five weighted subset identities agree')

    maxima = max_products_by_cardinality(NMAX)
    for n in range(1, NMAX + 1):
        m, s = triangular_witness(n)
        max_product = max(maxima[n])
        assert maxima[n][m] >= prod(s)
        for size, value in enumerate(maxima[n]):
            if value:
                assert size * (size + 1) <= 2 * n
                assert value * size ** size <= n ** size
        for r in range(1, 6):
            assert factorial(m) ** r <= prod(s) ** r <= arrays[r][n]
            assert arrays[r][n] <= arrays[0][n] * max_product ** r
    print(f'EXACT INTEGER BOUNDS: n=1..{NMAX}, r=1..5; triangular witness, AM-GM and count-times-max bounds agree')

    for n in range(2, NMAX + 1):
        root = sqrt_int(2 * n)
        log_root = times(Q(1, 2), log_int(2 * n))
        # q(n) <= exp(2 sqrt(n)).
        assert log_int(arrays[0][n])[1] <= 2 * sqrt_int(n)[0]
        # log Pmax <= R(log(n/R)+1) = R(log R-log 2+1).
        max_upper = mul(root, add(add(log_root, neg(LOG2)), const(1)))
        assert log_int(max(maxima[n]))[1] <= max_upper[0]
        m, _ = triangular_witness(n)
        if m > 1:
            stirling_lower = add(times(m, log_int(m)), const(1 - m))
            assert stirling_lower[1] <= log_int(factorial(m))[0]
        for r in range(1, 6):
            b = times(r, mul(root, add(log_root, const(-1))))
            lower = add(add(b, times(-2 * r, log_root)), const(-r))
            upper = add(add(b, times(r, mul(root, add(const(2), neg(LOG2))))), times(2, sqrt_int(n)))
            observed = log_int(arrays[r][n])
            assert lower[1] <= observed[0], ('lower', r, n)
            assert observed[1] <= upper[0], ('upper', r, n)
    print(f'CERTIFIED RATIONAL LOG BOUNDS: n=2..{NMAX}, r=1..5; counting, maximum-product, factorial and final logarithmic sandwiches hold')
    print(f'RATIONAL ENCLOSURES: atanh series {TERMS} terms; fixed denominator {SCALE}; square roots by integer isqrt')
    print('SCOPE: finite exact checks support the paper outline; no infinite limit or Lean proof is certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
