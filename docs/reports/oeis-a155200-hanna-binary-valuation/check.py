#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A155200 valuation conjecture.

Only Python's integer and Fraction arithmetic is used. Finite checks do not
prove the universal assertion. The independent exponential expansion checks
that the recurrence computes the source's ordinary coefficients.
"""

from fractions import Fraction
from math import factorial


LIMIT = 128
DIRECT_LIMIT = 24
SOURCE_PREFIX = (
    1, 2, 10, 188, 16774, 6745436, 11466849412, 80444398636280,
    2306003967992402758, 268654794629082985019564,
    126765597346260977505891041836,
    241678070948246232010898235031930952,
    1858395916567787793818891330877931472153500,
    57560683587056536617649234722821582390470430186648,
)


def recurrence_coefficients(limit):
    """Jovovic's recurrence in A155200, evaluated over the rationals."""
    coeffs = [Fraction(1)]
    for n in range(1, limit + 1):
        coeffs.append(sum((2 ** (k * k) * coeffs[n - k]
                           for k in range(1, n + 1)), Fraction(0)) / n)
    return coeffs


def multiply_truncated(left, right, limit):
    return [sum((left[k] * right[n - k] for k in range(n + 1)),
                Fraction(0)) for n in range(limit + 1)]


def exponential_coefficients(limit):
    """Directly expand sum_j H(x)^j/j!, without the coefficient recurrence."""
    exponent = [Fraction(0)] + [Fraction(2 ** (k * k), k)
                                for k in range(1, limit + 1)]
    power = [Fraction(1)] + [Fraction(0)] * limit
    result = list(power)
    for j in range(1, limit + 1):
        power = multiply_truncated(power, exponent, limit)
        result = [u + v / factorial(j) for u, v in zip(result, power)]
    return result


def normalized_coefficients(limit):
    """Independently evaluate the integral falling-factorial recurrence."""
    facts = [factorial(n) for n in range(limit + 1)]
    coeffs = [1]
    for n in range(1, limit + 1):
        coeffs.append(sum(facts[n - 1] // facts[n - k]
                          * 2 ** (k * k - k) * coeffs[n - k]
                          for k in range(1, n + 1)))
    return coeffs


def v2(positive_integer):
    assert positive_integer > 0
    return (positive_integer & -positive_integer).bit_length() - 1


def main():
    coeffs = recurrence_coefficients(LIMIT)
    assert coeffs[:len(SOURCE_PREFIX)] == list(SOURCE_PREFIX)
    print('FIDELITY: OEIS A155200 %S/%T/%U prefix n=0..13 matches exactly')
    assert exponential_coefficients(DIRECT_LIMIT) == coeffs[:DIRECT_LIMIT + 1]
    print(f'DIRECT_EXP: rational truncated exponential matches recurrence n=0..{DIRECT_LIMIT}')

    normalized = normalized_coefficients(LIMIT)
    for n, rational_value in enumerate(coeffs):
        assert rational_value.denominator == 1
        value = rational_value.numerator
        assert value > 0
        binary_sum = n.bit_count()
        assert v2(value) == binary_sum
        odd_quotient = Fraction(value, 2 ** binary_sum)
        assert odd_quotient.denominator == 1
        assert odd_quotient.numerator % 2 == 1
        b = rational_value * factorial(n) / 2 ** n
        assert b.denominator == 1
        assert b == normalized[n]
        assert normalized[n] % 4 == 1
        assert v2(factorial(n)) == n - binary_sum
        assert v2(value) + v2(factorial(n)) == n

    print(f'IDENTITY: positive integral a(n), v2(a(n))=s2(n), odd integral quotient n=0..{LIMIT}')
    print(f'NORMALIZATION: n! a(n)/2^n equals independent integral recurrence and is 1 mod 4 n=0..{LIMIT}')
    print(f'LEGENDRE: v2(n!)=n-s2(n) n=0..{LIMIT}')
    print('a(0..8) =', [int(x) for x in coeffs[:9]])
    print('v2(a(0..8)) =', [v2(int(x)) for x in coeffs[:9]])
    print('boundary (n, v2(a(n)), s2(n)) =',
          [(n, v2(int(coeffs[n])), n.bit_count())
           for n in (0, 1, 31, 32, 63, 64, 127, 128)])
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
