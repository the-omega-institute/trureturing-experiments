#!/usr/bin/env python3
"""Replay the numerical inequalities in the two RH continuation reports."""

from fractions import Fraction
from math import isqrt

SCALE = 10**24


def require(condition, label):
    if not condition:
        raise ValueError(label)


def ceiling(value):
    return -(-value.numerator // value.denominator)


def prime_sieve(limit):
    flags = bytearray(b"\1") * (limit + 1)
    flags[:2] = b"\0\0"
    for prime in range(2, isqrt(limit) + 1):
        if flags[prime]:
            flags[prime**2::prime] = bytes((limit - prime**2) // prime + 1)
    return flags


def integral_series(degree, denominator):
    total = term = Fraction(1)
    for order in range(degree):
        term *= Fraction(degree - order, denominator)
        total += term
    return total


def verify_quartic(primes):
    endpoint = 200_000
    require(3**11 < endpoint < Fraction(27, 10)**13, "quartic log brackets")
    epsilon = Fraction(13**2, 24 * isqrt(endpoint))
    terminal = (
        Fraction(5625, 2048) * Fraction(243, 241)**25
        * Fraction(endpoint, (endpoint - 1)**4)
        * ((1 + epsilon) * integral_series(25, 33) / 33 + 2 * epsilon / 11)
    )
    upper = ceiling(SCALE * terminal)
    require(upper == 44_446_336, "quartic terminal")
    cutoffs = {
        8231: (16952, 888884725255, Fraction(1, 2500000), 450000, Fraction(1, 10**12)),
        2543: (17612, 49669053146940, Fraction(1, 33750), 590500, Fraction(3, 10**7)),
        449: (17897, 24477013422083484, Fraction(7235955529, 450000000000),
              Fraction(82461631627375, 127401984), Fraction(1, 5000)),
    }
    choices = [Fraction(a, b) for a, b in ((1, 4), (3, 10), (1, 3), (3, 8), (2, 5))]
    count = 0
    for prime in range(endpoint, 448, -1):
        if prime in cutoffs:
            expected_count, expected_upper, mass, moment, margin = cutoffs[prime]
            require((count, upper) == (expected_count, expected_upper), f"quartic row {prime}")
            require(mass - moment * Fraction(upper, SCALE) > margin, f"quartic margin {prime}")
            print(f"Quartic cutoff {prime} gives live mass greater than {margin}")
        if primes[prime]:
            count += 1
            depth = Fraction(1, prime - 1)
            excess = 15*depth + 50*depth**2 + 60*depth**3 + 24*depth**4
            upper = min(ceiling(
                SCALE * 27 * depth**4 / (256 * delta**3 * (1 - delta))
                + upper * (1 + excess / (1 - delta))
            ) for delta in choices)


def quadratic_floor(value, prime):
    left = SCALE - value
    right = left + 4 * (3 * prime - 1) * value
    return value + (left + right + isqrt(4 * left * right)) // (4 * (prime - 1)**2)


def verify_quadratic(primes, capture=6997):
    endpoint = 5_000_000
    require(Fraction(11, 4)**15 < endpoint < Fraction(27, 10)**16, "quadratic log brackets")
    epsilon = Fraction(16**2, 24 * isqrt(endpoint))
    terminal = (
        Fraction(451, 449)**6 * Fraction(endpoint, (endpoint - 1)**2)
        * ((1 + epsilon) * integral_series(5, 15) + 2 * epsilon) / 15
    )
    lower, upper = 0, ceiling(SCALE * terminal)
    require(upper == 19_978_401_471_981_184, "quadratic terminal")
    count = 0
    captured = None
    for prime in range(endpoint, 73, -1):
        if prime == capture:
            captured = upper
        if prime == 6997:
            prefix_tail = upper
        if primes[prime]:
            count += 1
            lower = quadratic_floor(lower, prime)
            upper = quadratic_floor(upper, prime) + 1
    require(count == 348492, "quadratic prime count")
    require(lower == 7181063643367434956708, "quadratic lower replay")
    require(upper == 7187192385814385692968, "quadratic upper replay")
    require(Fraction(1391, 10) * upper < SCALE, "sufficient moment at 73")
    require(Fraction(6963, 50) * lower > SCALE, "scalar barrier at 73")
    require(upper < SCALE // 4, "quadratic monotonicity domain")
    require(prefix_tail == 27144070771294538066, "continuation at 6997")
    print("The quadratic continuation accepts 139.1 and has threshold below 139.26")
    require(captured is not None, "capture is outside the verified continuation")
    return captured


def main():
    from verify_minimum_modulus_nine import verify

    primes = prime_sieve(5_000_000)
    verify_quartic(primes)
    verify(verify_quadratic(primes))


if __name__ == "__main__":
    main()
