#!/usr/bin/env python3
"""Exact arithmetic checks for Spezia's OEIS A320043 square-value conjecture."""
from fractions import Fraction
from math import gcd, isqrt, prod
import sympy as sp


# OEIS A320043, revision 43: the entire DATA field, offset 1.
SOURCE_PREFIX = [
    1, 6, 13, 50, 37, 196, 189, 384, 351, 1210, 601, 2366, 1471,
    2156, 2941, 6936, 3277, 10830, 5563, 9022, 9681, 23276, 9897,
    26300, 19267, 30030, 23043, 58870, 21087, 76880, 46717, 59296,
    57801, 83546, 50281, 156066, 90973, 117968, 90539, 235340,
    86179, 284746,
]


def row(n):
    """Literal sum; check exact divisibility before each integer division."""
    assert n >= 1
    total = 0
    for k in range(1, n + 1):
        u = n + 1 - k
        d = gcd(u, k)
        assert u > 0 and k > 0 and d > 0
        numerator, denominator = u * u * k, d**3
        quotient, remainder = divmod(numerator, denominator)
        assert remainder == 0
        assert quotient == (u // d)**2 * (k // d)
        total += quotient
    return total


def rational_row(n):
    return sum((Fraction((n + 1 - k)**2 * k,
                         gcd(n + 1 - k, k)**3)
                for k in range(1, n + 1)), Fraction(0))


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    for n, expected in enumerate(SOURCE_PREFIX, 1):
        assert row(n) == expected, (n, row(n), expected)
    assert row(96) == 7376656 == 2716**2
    assert row(1) == 1 and row(6) == 196 == 14**2
    print(f'FIDELITY: all {len(SOURCE_PREFIX)} DATA terms agree with the literal sum')
    print('SOURCE_EXCEPTIONS: (1,1,1), (6,196,14), (96,7376656,2716)')

    # Independent rational evaluation checks the natural-division convention.
    for n in range(1, 257):
        assert rational_row(n) == row(n)
    small_squares = [(n, row(n), isqrt(row(n))) for n in range(1, 257)
                    if isqrt(row(n))**2 == row(n)]
    assert [n for n, _, _ in small_squares] == [1, 6, 96]
    print('SMALL_CASES: integer and rational row sums agree for every 1 <= n <= 256')
    print(f'SMALL_SQUARES: {small_squares}')

    p, k = sp.symbols('p k', integer=True)
    T = lambda j: (6*p**2*j*(j + 1) - 4*p*j*(j + 1)*(2*j + 1)
                   + 3*j**2*(j + 1)**2)
    assert sp.expand(T(k) - T(k - 1) - 12*k*(p - k)**2) == 0
    assert sp.expand(T(0)) == 0
    assert sp.expand(T(p - 1) - p**2*(p**2 - 1)) == 0
    print('TELESCOPE: exact symbolic difference and both endpoints verified')

    primes = [p for p in range(2, 258) if trial_prime(p)]
    for p in primes:
        assert all(gcd(p - k, k) == 1 for k in range(1, p))
        assert 12*row(p - 1) == p**2*(p**2 - 1)
    print(f'PRIME_ROWS: identity verified directly for all {len(primes)} primes p <= 257')

    x, y = sp.symbols('x y', integer=True)
    assert sp.expand((2*x*x - 1)**2 - 12*(2*x*y)**2 - 1
                     - 4*x*x*(x*x - 12*y*y - 1)) == 0
    chain = [(7, 2)]
    for _ in range(3):
        x, y = chain[-1]
        chain.append((2*x*x - 1, 2*x*y))
    assert chain == [(7, 2), (97, 28), (18817, 5432),
                     (708158977, 204427888)]
    assert all(x*x - 12*y*y == 1 for x, y in chain)
    assert 18817 == 31*607
    assert not trial_prime(18817)
    print(f'PELL_DOUBLING: {chain}; intermediate 18817 = 31 * 607 is composite')

    p, y = chain[-1]
    n, m = p - 1, p*y
    factors = {2: 9, 3: 1, 7: 2, 97: 2}
    assert prod(q**e for q, e in factors.items()) == p - 1
    assert all(trial_prime(q) for q in factors)
    assert pow(5, p - 1, p) == 1
    expected = [(2, 354079488, 708158976), (3, 236052992, 259204175),
                (7, 101165568, 342172932), (97, 7300608, 617955608)]
    for q, exponent, residue in expected:
        assert exponent == (p - 1)//q
        assert pow(5, exponent, p) == residue
        assert gcd(residue - 1, p) == 1
    print(f'LUCAS_CERTIFICATE: p-1 factors {factors}; base 5 residues {expected}')
    assert isqrt(p) == 26611
    assert trial_prime(p)
    print('INDEPENDENT_PRIMALITY: no integer divisor 2 through 26611')
    assert p*p == 501489136705686529
    assert 12*y*y == 501489136705686528
    assert n == 708158976 and n not in (1, 6, 96)
    assert m == 144767444036350576
    predicted = Fraction(p*p*(p*p - 1), 12)
    assert predicted.denominator == 1
    assert predicted.numerator == m*m == 20957612852817895952216852375531776
    print(f'WITNESS: n={n}, m={m}, prime-row value={predicted.numerator}=m^2')
    print('SCOPE: the large row is derived from the checked telescope and primality; not enumerated')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
