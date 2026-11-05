#!/usr/bin/env python3
"""Exact finite checks for OEIS A345096; these do not prove its limit."""
import argparse
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path


def literal(n):
    """The NAME formula, using exact rational ceiling/floor and integer powers."""
    total = Fraction(0)
    for k in range(1, n + 1):
        for i in range(1, k + 1):
            quotient = Fraction(n, k ** i)
            floor = quotient.numerator // quotient.denominator
            ceiling = -((-quotient.numerator) // quotient.denominator)
            exponent = 1 - ceiling + floor
            assert exponent in (0, 1)
            total += Fraction(k) ** exponent
    assert total.denominator == 1
    return total.numerator


def grouped(n):
    """Independent implementation of the proposed divisor-correction reduction."""
    excess = 0
    tau = sum(n % k == 0 for k in range(1, n + 1))
    for k in range(2, n + 1):
        t, power = 0, k
        while t < k and n % power == 0:
            t += 1
            power *= k
        if t:
            assert n % (k ** t) == 0
            assert 1 + t * (k - 1) <= k ** t <= n
        excess += (k - 1) * t
    assert 0 <= excess <= (n - 1) * (tau - 1)
    assert tau <= 2 * isqrt(n)
    value = n * (n + 1) // 2 + excess
    relative = Fraction(2 * value, n * n)
    residual = relative - 1 - Fraction(1, n)
    assert residual >= 0
    # Exactly checks residual <= 4/sqrt(n), with no irrational approximation.
    assert residual ** 2 <= Fraction(16, n)
    return value


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bfile', type=Path, help='OEIS A345096 b-file, n=1..1000')
    args = parser.parse_args()
    values = {n: grouped(n) for n in range(1, 1001)}
    assert [values[n] for n in range(1, 11)] == [1, 4, 8, 15, 19, 29, 34, 48, 57, 69]
    assert literal(4) == 15
    print('FIDELITY: source example a(4)=15 and first 10 entries agree')
    for n in range(1, 81):
        assert literal(n) == values[n], n
    print('IDENTITY: literal floor/ceiling sum equals grouped reduction for n=1..80')
    primes = [n for n in values if prime(n)]
    for p in primes:
        assert values[p] == (p * p + 3 * p - 2) // 2, p
    print(f'FIDELITY: source prime formula holds at all {len(primes)} primes <=1000')
    if args.bfile:
        table = {}
        raw = args.bfile.read_bytes()
        for line in raw.decode().splitlines():
            if line.strip() and not line.lstrip().startswith('#'):
                n, v = map(int, line.split())
                assert n not in table
                table[n] = v
        assert set(table) == set(range(1, 1001))
        assert table == values
        print('FIDELITY: all 1000 OEIS b-file entries agree')
        print('BFILE_SHA256=' + sha256(raw).hexdigest())
    print('BOUNDS: grouped excess, divisor pairing, and squared relative-error bound pass for n=1..1000')
    print('a(1000)=' + str(values[1000]))
    print('a(1000)/(1000^2/2)=' + str(Fraction(values[1000], 500000)))
    print('Finite checks only; the all-integer limit requires the stated general argument.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
