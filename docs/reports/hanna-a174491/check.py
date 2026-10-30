#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A174491 denominator conjecture."""

from fractions import Fraction
from math import factorial


def valuation_integer(n: int, p: int) -> int:
    assert n != 0 and p >= 2
    n = abs(n)
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def valuation(q: Fraction, p: int) -> int:
    assert q != 0
    return valuation_integer(q.numerator, p) - valuation_integer(q.denominator, p)


def recurrence_coefficients(limit: int) -> list[Fraction]:
    # x F' = F * Sum_{m>=1} x^m/2^(m^2), F(0) = 1.
    coefficients = [Fraction(1)]
    for n in range(1, limit + 1):
        coefficients.append(
            sum((coefficients[n - m] / (2 ** (m * m))
                 for m in range(1, n + 1)), Fraction(0)) / n
        )
    return coefficients


def product_coefficients(limit: int) -> list[Fraction]:
    # Independently multiply exp(x^m/(m*2^(m^2))) as truncated series.
    coefficients = [Fraction(1)] + [Fraction(0)] * limit
    for m in range(1, limit + 1):
        factor = [Fraction(0)] * (limit + 1)
        term = Fraction(1)
        factor[0] = term
        for b in range(1, limit // m + 1):
            term /= m * (2 ** (m * m)) * b
            factor[m * b] = term
        coefficients = [
            sum((coefficients[n - j] * factor[j]
                 for j in range(0, n + 1, m)), Fraction(0))
            for n in range(limit + 1)
        ]
    return coefficients


def partitions(n: int, least: int = 1):
    if n == 0:
        yield ()
    for m in range(least, n + 1):
        for suffix in partitions(n - m, m):
            yield (m,) + suffix


def partition_term(parts: tuple[int, ...]) -> Fraction:
    counts: dict[int, int] = {}
    for m in parts:
        counts[m] = counts.get(m, 0) + 1
    denominator = 1
    for m, b in counts.items():
        denominator *= (m * (2 ** (m * m))) ** b * factorial(b)
    return Fraction(1, denominator)


def main() -> None:
    coefficients = recurrence_coefficients(80)
    source_denominators = [
        1, 2, 32, 512, 262144, 33554432, 137438953472,
        562949953421312, 147573952589676412928,
        2417851639229258349412352, 2535301200456458802993406410752,
        2658455991569831745807614120560689152,
    ]
    source_numerators = [
        1, 1, 5, 19, 1921, 42253, 26628779, 15317124535,
        521786071318529, 1038675077390118457, 124715000994291451743437,
        14203783261714481442742242211, 49356730516809227634074385356860075,
    ]
    assert [c.denominator for c in coefficients[:12]] == source_denominators
    assert [c.numerator for c in coefficients[:13]] == source_numerators
    assert coefficients[:5] == [Fraction(1), Fraction(1, 2), Fraction(5, 32),
                                Fraction(19, 512), Fraction(1921, 262144)]
    print("FIDELITY: A174491 denominators n=0..11, A174490 numerators n=0..12, displayed series n=0..4: OK")

    product = product_coefficients(24)
    assert product == coefficients[:25]
    print("INDEPENDENT_EXPANSION: recurrence = product of exponential factors, n=0..24: OK")

    assert coefficients[0] == 1 and coefficients[0].denominator == 1
    for n in range(1, 81):
        assert coefficients[n] > 0
        assert coefficients[n].denominator == 2 ** (n * n + valuation_integer(n, 2))
        assert valuation(coefficients[n], 2) == -n * n - valuation_integer(n, 2)
    print("DENOMINATOR_IDENTITY: every n=1..80, plus constant term: OK")
    print("DENOMINATOR_EXPONENTS n=0..20:",
          [c.denominator.bit_length() - 1 for c in coefficients[:21]])

    partition_count = 0
    min_gap = None
    for n in range(1, 25):
        total = Fraction(0)
        distinguished = -n * n - valuation_integer(n, 2)
        for parts in partitions(n):
            partition_count += 1
            term = partition_term(parts)
            total += term
            if parts == (n,):
                assert valuation(term, 2) == distinguished
            else:
                distinct = len(set(parts))
                squares = sum(m * m for m in parts)
                extra = valuation_integer(term.denominator, 2) - squares
                assert squares <= n * n - 2 * n + 2
                assert extra <= n - distinct <= n - 1
                gap = valuation(term, 2) - distinguished
                assert gap > 0
                min_gap = gap if min_gap is None else min(min_gap, gap)
        assert total == coefficients[n]
    print(f"PARTITION_VALUATION: {partition_count} partitions for n=1..24, unique minimum at (n), minimum gap={min_gap}: OK")

    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    checks = 0
    for p in primes:
        for m in range(1, 129):
            h = Fraction(p, m * (2 ** (m * m)))
            if m % p == 0:
                k = m // p
                h -= Fraction(1, k * (2 ** (k * k)))
            assert h != 0 and valuation(h, p) >= 1
            checks += 1
        for k in range(1, 129):
            t = valuation_integer(k, p)
            # Check the required LTE lower bound by exact modular arithmetic.
            assert pow(2, (p * p - 1) * k * k, p ** (1 + 2 * t)) == 1
    print(f"ODD_PRIME_LOG_DIFFERENCE: {checks} coefficients, p=3..31 prime and m=1..128; LTE congruences k=1..128: OK")
    print("SCOPE: finite exact evidence; no proof of the unbounded conjecture is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
