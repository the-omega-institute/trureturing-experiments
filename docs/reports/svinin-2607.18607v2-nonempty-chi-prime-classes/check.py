#!/usr/bin/env python3
"""Exact convention and identity checks for Svinin arXiv:2607.18607v2.

This finite check supplies no proof of infinitude or density.
Requires Python 3 and SymPy. No floating-point arithmetic is used.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from math import gcd, prod
import sys
import sympy
from sympy import bernoulli, factorint, primerange

LIMIT = 2000
SOURCE_FIRST_TEN = {
    2: [3, 11, 23, 47, 59, 71, 83, 107, 131, 167],
    4: [5, 29, 53, 149, 173, 197, 269, 293, 317, 389],
    6: [7, 31, 67, 79, 103, 139, 151, 223, 283, 367],
    12: [13, 229, 277, 349, 373, 709, 733, 853, 877, 997],
    16: [17, 113, 593, 977, 1553, 2129, 2417, 2609, 2897, 3089],
    18: [19, 199, 307, 523, 739, 1063, 1171, 1279, 1423, 1531],
    36: [37, 397, 613, 829, 1117, 1549, 1693, 2557, 3637, 3709],
    40: [41, 281, 521, 761, 1481, 1721, 2441, 3881, 5081, 6521],
    42: [43, 211, 463, 547, 967, 1051, 1303, 1723, 2731, 3067],
    60: [61, 1021, 1381, 1741, 2221, 3181, 3541, 4021, 4261, 5821],
    72: [73, 1657, 2089, 2953, 5113, 5689, 7417, 8713, 9433, 10009],
    96: [97, 4129, 4513, 5857, 9697, 10273, 11617, 12577, 15073, 15649],
    100: [101, 701, 1301, 1901, 3701, 6101, 6701, 7901, 10301, 13901],
}
SOURCE_MINIMA_IN_RANGE = {8: 89, 20: 461, 22: 1871, 24: 409, 32: 353, 48: 1489}


@lru_cache(maxsize=None)
def B(n):
    # SymPy uses B_1=+1/2 in its one-argument API; fix that source convention.
    if n == 1:
        return Fraction(-1, 2)
    value = bernoulli(n)
    a, d = value.as_numer_denom()
    return Fraction(int(a), int(d))


def chi_direct(n):
    assert n >= 3 and n % 2 == 1
    b = B(n - 1)
    return Fraction((b / (n - 1)).denominator,
                    n * (b + Fraction(1, n)).denominator)


def F(m):
    assert m >= 2 and m % 2 == 0
    return prod(int(q) ** int(e) for q, e in factorint(m).items()
                if m % (int(q) - 1) == 0)


def valuation_integer(a, q):
    assert a != 0 and q >= 2
    a = abs(a)
    v = 0
    while a % q == 0:
        a //= q
        v += 1
    return v


def valuation_rational(a, q):
    return valuation_integer(a.numerator, q) - valuation_integer(a.denominator, q)


def main():
    assert sys.version_info >= (3, 8)
    print('Python=' + sys.version.split()[0] + ' SymPy=' + sympy.__version__)
    assert [B(n) for n in range(5)] == [Fraction(1), Fraction(-1, 2),
                                       Fraction(1, 6), Fraction(0), Fraction(-1, 30)]
    print('BERNOULLI_CONVENTION_OK B_0=1 B_1=-1/2 B_2=1/6 B_3=0 B_4=-1/30')
    primes = list(primerange(3, LIMIT))
    by_class = defaultdict(list)
    valuation_checks = 0
    for p in primes:
        m = p - 1
        b = B(m)
        direct = chi_direct(p)
        source_gcd = Fraction(m, gcd(m, abs(b.numerator)))
        product = F(m)
        assert direct == source_gcd == product, (p, direct, source_gcd, product)
        assert direct.denominator == 1 and direct.numerator > 0
        assert direct.numerator % 2 == 0
        by_class[direct.numerator].append(p)
        divided = b / m
        for q, e in factorint(m).items():
            q, e = int(q), int(e)
            v = valuation_rational(divided, q)
            if m % (q - 1) == 0:
                assert v == -(e + 1), (p, q, e, v)
            else:
                assert v >= 0, (p, q, e, v)
            valuation_checks += 1
        print('p=' + str(p) + ' chi_direct=' + str(direct) +
              ' source_gcd=' + str(source_gcd) + ' product=' + str(product))
    source_table_checks = 0
    complete_rows = 0
    for mu, source_row in SOURCE_FIRST_TEN.items():
        expected = [p for p in source_row if p < LIMIT]
        actual = by_class[mu][:len(expected)]
        assert actual == expected, (mu, actual, expected)
        if len(expected) < 10:
            assert by_class[mu] == expected, (mu, by_class[mu], expected)
        else:
            complete_rows += 1
        source_table_checks += len(expected)
        print('SOURCE_TABLE_OK mu=' + str(mu) + ' prefix=' + ','.join(map(str, actual)))
    for mu, minimum in SOURCE_MINIMA_IN_RANGE.items():
        assert by_class[mu][0] == minimum, (mu, minimum, by_class[mu][0])
        print('SOURCE_MINIMUM_OK mu=' + str(mu) + ' minimum=' + str(minimum))
    print('DOMAIN all odd primes 3<=p<' + str(LIMIT) + '; p=2 is outside the source domain')
    print('IDENTITY_CHECKS=' + str(len(primes)) + ' MAX_PRIME=' + str(primes[-1]))
    print('DIVIDED_BERNOULLI_VALUATION_CHECKS=' + str(valuation_checks))
    print('SOURCE_TABLE_ENTRIES=' + str(source_table_checks) + ' COMPLETE_ROWS=' + str(complete_rows))
    print('SOURCE_MINIMA=' + str(len(SOURCE_MINIMA_IN_RANGE)))
    print('FINITE_CHECK_ONLY infinitude and density remain unproved by this program')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
