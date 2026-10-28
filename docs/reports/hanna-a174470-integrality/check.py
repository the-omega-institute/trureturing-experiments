#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A174470 integrality conjecture."""

from fractions import Fraction
from hashlib import sha256
from math import factorial

LIMIT = 120
# OEIS A174470, %S/%T/%U, indices 0..26, entry revision #2.
SOURCE_PREFIX = [
    1, 3, 9, 18, 54, 162, 405, 1215, 3645, 11907, 35721, 107163,
    279936, 839808, 2519424, 6948099, 20844297, 62532891, 199762767,
    599288301, 1797864903, 4931772480, 14795317440, 44385952320,
    125801650638, 377404951914, 1132214855742,
]


def valuation(p, n):
    assert n > 0 and p >= 2
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e


def b(n):
    return int(n > 0 and valuation(3, n) % 2 == 0)


def exp_coefficients(base, limit):
    """n*a[n] = sum_{j=1}^n b(j)*base^j*a[n-j], over Q."""
    a = [Fraction(1)]
    for n in range(1, limit + 1):
        a.append(sum((b(j) * base**j * a[n-j] for j in range(1, n+1)),
                     Fraction(0)) / n)
    return a


def multiply(a, f, limit):
    out = [Fraction(0) for _ in range(limit + 1)]
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(f[:limit+1-i]):
                if y:
                    out[i+j] += x*y
    return out


def fractional_euler_product(limit):
    out = [Fraction(1)] + [Fraction(0)] * limit
    d, k, factors = 1, 0, 0
    while d <= limit:
        s = (-1)**(k+1)
        alpha = Fraction(s, d)
        factor = [Fraction(0)] * (limit+1)
        factor[0] = Fraction(1)
        binomial = Fraction(1)
        for m in range(1, limit//d+1):
            binomial *= (alpha - (m-1)) / m
            coefficient = binomial * (-3**d)**m
            assert coefficient.denominator == 1, (k, m, coefficient)
            numerator_product = 1
            for j in range(m):
                numerator_product *= s-d*j
            explicit = Fraction((-1)**m * 3**((d-k)*m)
                                * numerator_product, factorial(m))
            assert coefficient == explicit
            factor[d*m] = coefficient
        out = multiply(out, factor, limit)
        factors += 1
        k += 1
        d *= 3
    return out, factors


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, p) if d*d <= p)


def main():
    # Independent source definition: the morphism 1->110, 0->111.
    word = [1]
    for _ in range(7):
        word = [digit for x in word for digit in ([1, 1, 0] if x else [1, 1, 1])]
    assert word == [b(n) for n in range(1, len(word)+1)]
    for n in range(1, 10001):
        assert b(3*n) == 1-b(n)
        assert b(3*n+1) == b(3*n+2) == 1
        alternating = 0
        d, k = 1, 0
        while d <= n:
            alternating += (-1)**k * (n//d - (n-1)//d)
            d *= 3
            k += 1
        assert b(n) == alternating
    print(f'A014578 fidelity: morphism matches indices 1..{len(word)}; '
          'recurrence and floor/divisibility formula match n=1..10000')

    a = exp_coefficients(3, LIMIT)
    assert a[:len(SOURCE_PREFIX)] == SOURCE_PREFIX
    assert all(x.denominator == 1 for x in a)
    print('A174470 fidelity: all 27 displayed coefficients match, indices 0..26')
    print('source_prefix=' + ','.join(map(str, SOURCE_PREFIX)))
    print(f'Integrality check: every rational coefficient at indices 0..{LIMIT} '
          'has denominator 1')
    product, factor_count = fractional_euler_product(LIMIT)
    assert product == a
    print(f'Fractional Euler product: {factor_count} locally relevant factors '
          f'match exponential coefficients at indices 0..{LIMIT}')

    primes = [p for p in range(2, 102) if is_prime(p)]
    for p in primes:
        for m in range(1, 301):
            modulus = p**(valuation(p, m)+1)
            difference = b(p*m)*pow(3, p*m, modulus)-b(m)*pow(3, m, modulus)
            assert difference % modulus == 0, (p, m)
    print(f'Dold congruence check: {len(primes)} primes p<=101 and 1<=m<=300; '
          f'{len(primes)*300} cases')

    for base in range(-12, 13):
        vals = exp_coefficients(base, 45)
        assert vals[3] == Fraction(2*base**3, 3)
        assert all(x.denominator == 1 for x in vals) == (base % 3 == 0)
        if base % 3 == 0:
            assert vals == [a[n]*(base//3)**n for n in range(46)]
    print('Sharp-extension finite check: B=-12..12, coefficients 0..45; '
          'integral iff 3 divides B; [x^3]=2*B^3/3')
    payload = '\n'.join(str(x) for x in a).encode('ascii')
    print('coefficients_0_120_sha256=' + sha256(payload).hexdigest())
    print('Finite checks are evidence for the planned argument, not an all-index proof.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
