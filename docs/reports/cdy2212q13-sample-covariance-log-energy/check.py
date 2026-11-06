#!/usr/bin/env python3
"""Exact checks for the sample-covariance logarithmic-energy preregistration.

The small-noise coefficient L_n is the coefficient of log(1/e) in the
pairwise energy for the Gaussian-regularised Bernoulli law.  For a Bernoulli
matrix A and an independent copy C, the determinant identity gives

    k(A,C) = 2 z(A) z(C) + sum_lambda m_A(lambda) m_C(lambda),

where z is the nullity of A and the sum is over common positive Gram
eigenvalues.  The order-one contribution is generic for nonzero clusters;
the zero clusters have width e^2.  All quantities below are exact integers,
SymPy rationals, or fractions.Fraction values.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product

import sympy as sp


X = sp.symbols("x")


def source_laguerre_formula(n: sp.Integer | int) -> sp.Expr:
    """The source Theorem 1.3 formula for the square Gaussian case."""
    n = sp.sympify(n)
    return sp.Rational(3, 2) + sp.log(n) + sp.EulerGamma + sp.Rational(1, 2) / n - sp.harmonic(n)


def matrix_signature(n: int, bits: tuple[int, ...]) -> tuple:
    """Return exact nullity and positive Gram characteristic factors.

    The two Bernoulli values are represented after the common nonzero scale
    by -1 and 255, so equality of Gram eigenvalues is unchanged.
    """
    entries = [-1 if bit == 0 else 255 for bit in bits]
    matrix = sp.Matrix(n, n, entries)
    gram = matrix * matrix.T
    factors = sp.factor_list(gram.charpoly(X).as_expr())[1]
    nullity = 0
    positive_factors: list[tuple[tuple[int, ...], int]] = []
    for factor, multiplicity in factors:
        polynomial = sp.Poly(factor, X, domain=sp.QQ).monic()
        if polynomial == sp.Poly(X, X, domain=sp.QQ):
            nullity += multiplicity
        else:
            positive_factors.append((tuple(int(c) for c in polynomial.all_coeffs()), multiplicity))
    return (sum(bits), nullity, tuple(sorted(positive_factors)))


def signatures(n: int) -> Counter:
    return Counter(
        matrix_signature(n, bits)
        for bits in product((0, 1), repeat=n * n)
    )


def order(signature_a: tuple, signature_c: tuple) -> int:
    _, nullity_a, factors_a = signature_a
    _, nullity_c, factors_c = signature_c
    common_nonzero = sum(
        multiplicity_a * multiplicity_c
        for factor_a, multiplicity_a in factors_a
        for factor_c, multiplicity_c in factors_c
        if factor_a == factor_c
    )
    return 2 * nullity_a * nullity_c + common_nonzero


def exact_L(n: int, p: Fraction) -> tuple[Fraction, Counter]:
    """Compute L_n exactly by enumerating the finite Bernoulli support."""
    grouped = signatures(n)
    q = 1 - p
    total = Fraction(0)
    order_mass: Counter = Counter()
    for signature_a, count_a in grouped.items():
        for signature_c, count_c in grouped.items():
            count_b = signature_a[0] + signature_c[0]
            weight = (
                count_a
                * count_c
                * p**count_b
                * q ** (2 * n * n - count_b)
            )
            k = order(signature_a, signature_c)
            total += weight * k
            order_mass[k] += weight
    return total / (n * n), order_mass


def main() -> None:
    # Fidelity: the source's square-Laguerre formula has MP limit 3/2.
    assert source_laguerre_formula(1) == 1 + sp.EulerGamma
    limit_value = sp.simplify(sp.limit(source_laguerre_formula(sp.Symbol("n", positive=True)), sp.Symbol("n", positive=True), sp.oo))
    assert limit_value == sp.Rational(3, 2)
    print("source Gaussian formula: n=1 -> 1 + EulerGamma")
    print(f"source Gaussian formula limit: {limit_value}")

    p = Fraction(1, 256)
    l2, mass2 = exact_L(2, p)
    l3, mass3 = exact_L(3, p)
    all_common_2 = matrix_signature(2, (0,) * 4)
    all_common_3 = matrix_signature(3, (0,) * 9)
    assert order(all_common_2, all_common_2) == 3
    assert order(all_common_3, all_common_3) == 9
    assert l2 <= Fraction(101, 128)
    assert l3 >= Fraction(119, 128)
    assert l3 > l2
    print(f"L_2 = {l2}")
    print(f"L_3 = {l3}")
    print(f"L_3 - L_2 = {l3 - l2}")
    print(f"L_2 order masses = {sorted(mass2.items())}")
    print(f"L_3 order masses = {sorted(mass3.items())}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
