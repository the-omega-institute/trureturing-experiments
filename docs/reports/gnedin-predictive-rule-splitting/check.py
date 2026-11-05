#!/usr/bin/env python3
"""Exact checks for the all-degree predictive-rule splitting candidate.

The checks use SymPy expressions over QQ and Python Fractions only.  They
reproduce the signed cubic factorisation displayed in Gnedin, v2, and test
the triangular-zero and largest-moved-index calculations for every
permutation through d = 9 (46,233 permutations in total).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from math import factorial

import sympy as sp


U, V, X = sp.symbols("u v x")


def rising(value: Fraction | sp.Expr, length: int):
    result = 1
    for i in range(length):
        result *= value + i
    return result


def bernstein_coefficients(poly: sp.Expr, degree: int) -> list[sp.Rational]:
    """Return c_j with poly = sum c_j*x^j*(1-x)^(degree-j)."""
    cs = sp.symbols(f"c0:{degree + 1}")
    equation = sp.Poly(
        sp.expand(sum(cs[j] * X**j * (1 - X) ** (degree - j)
                     for j in range(degree + 1)) - poly), X
    )
    solution = sp.solve(equation.all_coeffs(), cs, dict=True)
    assert len(solution) == 1
    return [sp.Rational(solution[0][c]) for c in cs]


def q_from_coefficients(cs: list[sp.Rational], degree: int,
                        u: sp.Expr = U, v: sp.Expr = V) -> sp.Expr:
    return sp.expand(sum(cs[j] * rising(u, j) * rising(v, degree - j)
                         for j in range(degree + 1)))


def triangular_zero_checks() -> None:
    # The source's Lemma 6.13 consequence: Q(-r,-s)=0 whenever r+s<d.
    for degree in range(1, 13):
        cs = [sp.Rational(i + 1) for i in range(degree + 1)]
        for r in range(degree):
            for s in range(degree - r):
                value = q_from_coefficients(cs, degree, sp.Integer(-r), sp.Integer(-s))
                assert value == 0, (degree, r, s, value)


def cubic_factorisation_check() -> None:
    # The exceptional cofactor in Theorem 6.14 / Remark 6.15.
    cubic = sp.expand((2 - X) * (1 + X) * (1 - 2 * X))
    assert cubic == 2 - 3 * X - 3 * X**2 + 2 * X**3
    cs = bernstein_coefficients(cubic, 3)
    assert cs == [sp.Integer(2), sp.Integer(3), sp.Integer(-3), sp.Integer(-2)]
    q = q_from_coefficients(cs, 3)
    expected = (V - U) * (U + 2 * V + 2) * (2 * U + V + 2)
    assert sp.expand(q - expected) == 0
    # The signed cubic changes sign at x=1/2, so it is not a probability
    # density after multiplication by any nonzero scalar.
    assert cubic.subs(X, sp.Rational(1, 4)) > 0
    assert cubic.subs(X, sp.Rational(3, 4)) < 0


def largest_moved_index_checks() -> None:
    checked = 0
    for degree in range(2, 10):
        domain = tuple(range(1, degree))
        for perm in permutations(domain):
            checked += 1
            sigma = dict(zip(domain, perm))
            moved = [j for j in domain if sigma[j] != j]
            if not moved:
                continue
            J = max(moved)
            assert all(sigma[j] == j for j in domain if j > J)
            assert sigma[J] < J

            # At (-1, -(J-1)), every product factor is nonzero.  Work in
            # Fractions so this remains an exact sign check.
            for j in domain:
                value = Fraction(-j) + Fraction(-sigma[j] * (J - 1)) + j * sigma[j]
                assert value != 0, (degree, perm, J, j, value)

            # Reindexing by sigma and swapping u,v is the inverse-permutation
            # version of the same check at (-(J-1), -1).
            inverse = {value: key for key, value in sigma.items()}
            for j in domain:
                value = Fraction(-j * (J - 1)) + Fraction(-inverse[j]) + j * inverse[j]
                assert value != 0, (degree, perm, J, j, value)

            # The two grid equations force J = 2.  For J != 2 the linear
            # system has only alpha = beta = 0, which is forbidden.
            if J != 2:
                determinant = 1 - (J - 1) ** 2
                assert determinant != 0, (degree, perm, J)

    assert checked == 46233, checked


def homogeneous_shape_checks() -> None:
    # Identity permutations produce degree-at-most-one cofactors; the only
    # nonidentity survivor is the swap (1 2), whose cofactor is cubic.
    for degree in range(3, 10):
        identity = list(range(1, degree))
        factors = [j * U + identity[j - 1] * V for j in identity]
        h_identity = sp.expand((U + V) ** (degree - 1))
        assert sp.expand(sp.prod(factors) - factorial(degree - 1) * h_identity) == 0

        swap = [2, 1] + list(range(3, degree))
        h_swap = sp.expand((V - U) * (U + 2 * V) * (2 * U + V) * (U + V) ** (degree - 3))
        raw_swap = sp.expand((V - U) * sp.prod(
            j * U + swap[j - 1] * V for j in range(1, degree)
        ))
        # The omitted positive scalar is irrelevant to proportionality.
        assert sp.expand(raw_swap - (factorial(degree - 1) // 2) * h_swap) == 0
        p_swap = sp.expand(h_swap.subs({U: X, V: 1 - X}))
        assert sp.Poly(p_swap, X).degree() == 3


def main() -> None:
    triangular_zero_checks()
    cubic_factorisation_check()
    largest_moved_index_checks()
    homogeneous_shape_checks()
    print("triangular zeros: OK (degrees 1..12)")
    print("signed cubic factorisation and sign change: OK")
    print("largest-moved-index checks: OK (46,233 permutations, d=2..9)")
    print("homogeneous degree comparison: OK (d=3..9)")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
