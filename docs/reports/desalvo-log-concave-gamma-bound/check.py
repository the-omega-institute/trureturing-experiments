#!/usr/bin/env python3
"""Exact checks for the DeSalvo Open Problem 8.2 Gamma clause.

All evidence is Fraction arithmetic.  The last section checks the ratio
identity used in the source's Bessel Stein operator as a fidelity check.
"""

from fractions import Fraction
from math import gcd
import sys


def pmf(M: int, n: int) -> Fraction:
    """The proposed full-support log-concave family p_M(n)."""
    D = M * M + M - 1
    if n == 0:
        return Fraction(M - 1, D)
    return Fraction(M - 1, D) * Fraction(M) ** (2 - n)


def death(M: int, n: int) -> Fraction:
    if n == 0:
        return Fraction(0)
    return pmf(M, n - 1) / pmf(M, n)


def exact_moments(M: int):
    D = M * M + M - 1
    mean = Fraction(M**3, (M - 1) * D)
    second = Fraction(M**3 * (M + 1), (M - 1) ** 2 * D)
    variance = second - mean * mean
    stated_variance = Fraction(M**3 * (2 * M * M - 1), (M - 1) ** 2 * D * D)
    assert variance == stated_variance
    return mean, second, variance


def check_family(M: int) -> None:
    assert M >= 2
    D = M * M + M - 1

    # Normalization and positivity, with a finite head plus an exact geometric tail.
    p0 = pmf(M, 0)
    p1 = pmf(M, 1)
    q = Fraction(1, M)
    assert p0 > 0 and p1 > 0
    assert p0 + p1 / (1 - q) == 1
    assert all(pmf(M, n) > 0 for n in range(0, 12))

    # Consecutive ratios are M, 1/M, 1/M, ...; hence the pmf is log-concave.
    assert pmf(M, 1) / pmf(M, 0) == M
    assert pmf(M, 2) / pmf(M, 1) == q
    assert all(pmf(M, n + 1) / pmf(M, n) == q for n in range(2, 12))
    assert all(pmf(M, n) ** 2 >= pmf(M, n - 1) * pmf(M, n + 1) for n in range(1, 11))

    # The induced death rates have the two values claimed in the nomination.
    assert death(M, 1) == Fraction(1, M)
    assert all(death(M, n) == M for n in range(2, 12))

    # Directly verify the indicator Stein solution for A={0}, w=1.
    g = {0: Fraction(0), 1: Fraction(M * M, D)}
    for n in range(2, 14):
        g[n] = Fraction(1, D)
    pA = p0
    assert g[1] - death(M, 0) * g[0] == 1 - pA
    assert g[2] - death(M, 1) * g[1] == -pA
    for n in range(2, 13):
        assert g[n + 1] - death(M, n) * g[n] == -pA

    mean, second, variance = exact_moments(M)
    assert variance > 0
    gamma_lower = g[1]
    ratio_sq = gamma_lower * gamma_lower / variance
    claimed_lower = Fraction(M * (M - 1) ** 2, 2 * M * M - 1)
    assert ratio_sq == claimed_lower
    assert ratio_sq >= Fraction(M, 8)

    # The concrete values quoted in the nomination are exact at M=100.
    if M == 100:
        assert D == 10099
        assert p0 == Fraction(99, 10099)
        assert p1 == Fraction(9900, 10099)
        assert mean == Fraction(1000000, 999801)
        assert variance == Fraction(19999000000, 999602039601)
        assert ratio_sq == Fraction(980100, 19999)


def source_bessel_ratio_fidelity() -> None:
    """Exact fidelity check of the source's displayed Bessel ratio identity."""
    # The source defines unnormalised weights q_n = nu^n/(n!(n+k)!).
    k, nu = 3, Fraction(5)
    q = lambda n: nu**n / (factorial(n) * factorial(n + k))
    for n in range(0, 9):
        assert q(n + 1) / q(n) == nu / ((n + 1) * (n + k + 1))


def factorial(n: int) -> int:
    out = 1
    for i in range(2, n + 1):
        out *= i
    return out


def main() -> int:
    for M in (2, 3, 10, 100, 1000):
        check_family(M)
        print(f"M={M}: exact normalization, log-concavity, Stein solution, moments, and Gamma/sigma bound OK")
    source_bessel_ratio_fidelity()
    print("fidelity: source Bessel ratio identity checked exactly for k=3, nu=5")
    # For every finite C, choose M > max(2, 8 C^2); the exact ratio lower bound
    # then gives Gamma/sigma > C.  This is the quantified refutation mechanism.
    print("general bound: (Gamma/sigma)^2 >= M/8 for every integer M >= 2")
    print("ALL_CHECKS_OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
