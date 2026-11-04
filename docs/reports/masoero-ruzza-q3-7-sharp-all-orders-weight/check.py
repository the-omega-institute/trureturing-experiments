#!/usr/bin/env python3
"""Exact checks for Masoero--Ruzza Question 3.7.

All arithmetic is over SymPy's exact rational polynomial domains.  The Wronskian
is evaluated after removing its common Gaussian and monomial factors, using the
source's monic generalized Laguerre coefficients.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import permutations

import sympy as sp

T = sp.symbols("T")
Z = sp.symbols("Z")
B = sp.symbols("beta")


def conjugate(lam: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(part >= j for part in lam) for j in range(1, lam[0] + 1)) if lam else ()


def partitions(n: int, largest: int | None = None):
    if n == 0:
        yield ()
        return
    largest = n if largest is None else min(n, largest)
    for first in range(largest, 0, -1):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def laguerre_monic(n: int, alpha: sp.Rational) -> sp.Poly:
    # [T^i] \hat L_n^(alpha) = (-1)^(n+i) binom(n,i) (alpha+i+1)_(n-i).
    expr = sum(
        (-1) ** (n + i)
        * sp.binomial(n, i)
        * sp.rf(alpha + i + 1, n - i)
        * T**i
        for i in range(n + 1)
    )
    return sp.Poly(expr, T, domain=sp.QQ)


def vandermonde(vector: tuple[sp.Expr, ...]) -> sp.Expr:
    value = sp.Integer(1)
    for i in range(len(vector)):
        for j in range(i + 1, len(vector)):
            value *= vector[j] - vector[i]
    return sp.expand(value)


def permutation_determinant(matrix: list[list[sp.Poly]]) -> sp.Poly:
    size = len(matrix)
    result = sp.Poly(0, T, domain=sp.QQ)
    for perm in permutations(range(size)):
        inversions = sum(perm[i] > perm[j] for i in range(size) for j in range(i + 1, size))
        term = sp.Poly(1, T, domain=sp.QQ)
        for row, col in enumerate(perm):
            term *= matrix[row][col]
        result += (-1) ** inversions * term
    return result


@lru_cache(maxsize=None)
def phi_at(lam: tuple[int, ...], beta: sp.Rational) -> sp.Poly:
    """Return Phi_{lambda,lambda'}^(beta)(T), exactly, from the Wronskian."""
    mu = conjugate(lam)
    r, s = len(lam), len(mu)
    size = r + s
    shifted_beta = beta - r + s
    m = tuple(lam[i] + r - i - 1 for i in range(r))
    n = tuple(mu[i] + s - i - 1 for i in range(s))
    if size == 0:
        return sp.Poly(1, T, domain=sp.QQ)

    columns = [laguerre_monic(q, shifted_beta) for q in m]
    columns += [laguerre_monic(q, -shifted_beta) for q in n]
    exponents = [sp.Rational(1, 2) + shifted_beta] * r
    exponents += [sp.Rational(1, 2) - shifted_beta] * s

    # If f=e^(-x^2/2)x^a p(T), then D(f)=e^(-x^2/2)x^(a-1)
    # ((a-T)p+2T p').  This produces the Wronskian after its common
    # Gaussian and monomial factors have been divided out.
    rows: list[list[sp.Poly]] = []
    for derivative in range(size):
        rows.append(columns)
        columns = [
            sp.Poly(
                (a - derivative) * p.as_expr()
                + 2 * T * sp.diff(p.as_expr(), T)
                - T * p.as_expr(),
                T,
                domain=sp.QQ,
            )
            for p, a in zip(columns, exponents)
        ]

    determinant = permutation_determinant(rows)
    extra_power = (r - s) ** 2 + size**2 - 2 * size
    assert extra_power % 4 == 0
    extra_power //= 4
    if extra_power:
        determinant = determinant.exquo(sp.Poly(T**extra_power, T, domain=sp.QQ))

    kappa = 2 ** (size * (size - 1) // 2) * vandermonde(m) * vandermonde(n)
    for mi in m:
        for nj in n:
            kappa *= nj - mi - shifted_beta
    assert kappa != 0, (lam, beta, m, n, shifted_beta)
    result = determinant.exquo(sp.Poly(kappa, T, domain=sp.QQ))
    for (degree,), _ in result.terms():
        assert degree % 2 == 0, (lam, beta, result)
    return result


@lru_cache(maxsize=None)
def polynomial_at(lam: tuple[int, ...], beta: sp.Rational) -> sp.Poly:
    phi = phi_at(lam, beta)
    result = sp.Poly(0, Z, domain=sp.QQ)
    for (degree,), coefficient in phi.terms():
        result += sp.Poly(coefficient * Z ** (degree // 2), Z, domain=sp.QQ)
    assert result.monic()
    return result


def root_power_sums(poly: sp.Poly, maximum: int) -> list[sp.Rational]:
    degree = poly.degree()
    coefficients = poly.all_coeffs()
    powers: list[sp.Rational] = []
    for k in range(1, maximum + 1):
        value = sp.Integer(0)
        for j in range(1, min(k - 1, degree) + 1):
            value += coefficients[j] * powers[k - j - 1]
        if k <= degree:
            value += k * coefficients[k]
        powers.append(sp.factor(-value))
    return powers


def shifted_p(lam: tuple[int, ...], k: int) -> sp.Rational:
    if k == 0:
        return sp.Integer(1)
    return sp.expand(
        sum(
            (sp.Rational(part - i, 1) + sp.Rational(1, 2)) ** (k - 1)
            - (sp.Rational(-i, 1) + sp.Rational(1, 2)) ** (k - 1)
            for i, part in enumerate(lam, start=1)
        )
    )


def source_s1(lam: tuple[int, ...], beta: sp.Rational) -> sp.Expr:
    p2, p3, p4 = (shifted_p(lam, j) for j in (2, 3, 4))
    return sp.expand(2 * p4 - 2 * p2**2 + 3 * beta * p3 + beta**2 * p2 + p2 / 2)


def source_s2(lam: tuple[int, ...], beta: sp.Rational) -> sp.Expr:
    p2, p3, p4, p5, p6 = (shifted_p(lam, j) for j in (2, 3, 4, 5, 6))
    return sp.expand(
        14 * p6 - 40 * p4 * p2 + sp.Rational(64, 3) * p2**3
        + sp.Rational(65, 3) * p4 - 20 * p2**2 + sp.Rational(27, 8) * p2
        + beta * (35 * p5 - 60 * p3 * p2 + sp.Rational(65, 2) * p3)
        + beta**2 * (30 * p4 - 22 * p2**2 + sp.Rational(25, 2) * p2)
        + 10 * beta**3 * p3 + beta**4 * p2
    )


def interpolated_power(lam: tuple[int, ...], k: int, samples: list[sp.Rational]) -> sp.Poly:
    values = [(beta, root_power_sums(polynomial_at(lam, beta), k)[k - 1]) for beta in samples]
    return sp.Poly(sp.interpolate(values, B), B, domain=sp.QQ)


def main() -> None:
    small_partitions = [lam for n in range(5) for lam in partitions(n)]
    sample_betas = [sp.Rational(2 * j + 1, 2) for j in range(-8, 9)]

    # Fidelity with the source's stated S_1 and S_2 formulas, at exact rational beta.
    for lam in small_partitions:
        for beta in sample_betas:
            powers = root_power_sums(polynomial_at(lam, beta), 2)
            assert powers[0] == source_s1(lam, beta), ("S1", lam, beta, powers[0], source_s1(lam, beta))
            assert powers[1] == source_s2(lam, beta), ("S2", lam, beta, powers[1], source_s2(lam, beta))

    # Exact interpolation checks for k <= 4.  The interpolated polynomial is
    # formed from 17 rational samples, then evaluated at fresh samples.
    for lam in small_partitions:
        for k in range(1, 5):
            polynomial = interpolated_power(lam, k, sample_betas)
            assert polynomial.degree() <= 2 * k + 2, ("degree", lam, k, polynomial)
            for beta in (sp.Rational(19, 2), sp.Rational(-21, 2), sp.Rational(23, 2)):
                expected = polynomial.eval(beta)
                actual = root_power_sums(polynomial_at(lam, beta), k)[k - 1]
                assert actual == expected, ("continuation", lam, k, beta, actual, expected)

    print("Exact Wronskian fidelity: S_1 and the source S_2 formula for 15 partitions, 17 rational beta values.")
    print("Sharp filtration check: k=1,2,3,4 for the same partitions, with fresh exact beta evaluations.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
