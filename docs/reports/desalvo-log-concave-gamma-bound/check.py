#!/usr/bin/env python3
"""Exact checks for the DeSalvo Open Problem 8.2 Gamma clause.

All evidence is Fraction arithmetic.  The last section checks the ratio
identity used in the source's Bessel Stein operator as a fidelity check.
"""

from fractions import Fraction
import sys

import sympy as sp


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
    # Differentiate the exact probability generating function, including its
    # infinite geometric tail, rather than checking the proposed formulas alone.
    z = sp.Symbol("z")
    P = sp.Rational(M - 1, D) + sp.Rational(M * (M - 1), D) * z / (1 - z / M)
    mean = Fraction(sp.diff(P, z).subs(z, 1))
    second = Fraction((sp.diff(P, z, 2) + sp.diff(P, z)).subs(z, 1))
    assert mean == Fraction(M**3, (M - 1) * D)
    assert second == Fraction(M**3 * (M + 1), (M - 1) ** 2 * D)
    variance = second - mean * mean
    assert variance == Fraction(M**3 * (2 * M * M - 1), (M - 1) ** 2 * D * D)
    return mean, second, variance


def symbolic_family_identities() -> None:
    # A symbolic computation independent of the finite parameter checks.
    M, z = sp.symbols("M z", positive=True)
    D = M**2 + M - 1
    P = (M - 1) / D + M * (M - 1) * z / (D * (1 - z / M))
    mean = sp.diff(P, z).subs(z, 1)
    second = (sp.diff(P, z, 2) + sp.diff(P, z)).subs(z, 1)
    variance = sp.factor(second - mean**2)
    assert sp.factor(P.subs(z, 1) - 1) == 0
    assert sp.factor(mean - M**3 / ((M - 1) * D)) == 0
    assert sp.factor(variance - M**3 * (2 * M**2 - 1) / ((M - 1)**2 * D**2)) == 0
    ratio_sq = sp.factor((M**2 / D)**2 / variance)
    assert sp.factor(ratio_sq - M * (M - 1)**2 / (2 * M**2 - 1)) == 0
    # Clearing the positive denominator gives a polynomial with nonnegative
    # coefficients after M=t+2. This verifies the inequality for every M>=2.
    t = sp.Symbol("t", nonnegative=True)
    numerator = sp.factor((ratio_sq - M / 8) * 8 * (2 * M**2 - 1) / M)
    poly = sp.Poly(sp.expand(numerator.subs(M, t + 2)), t)
    assert poly.as_expr() == 6 * t**2 + 8 * t + 1
    assert all(c > 0 for c in poly.all_coeffs())
    print("symbolic family: exact PGF moments and ratio^2 >= M/8 for every M>=2")


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
    symbolic_family_identities()
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
