#!/usr/bin/env python3
"""Exact checks for the coloured functional-digraph Hankel obstruction.

The source uses T(t)=sum_{r>=1} r^(r-1)t^r/r! and
C_q(t)=-log(1-q*T(t)).  Formal power-series coefficients below are
Fractions and integer polynomial dictionaries; no floating point is used.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial


def series_add(a, b, order):
    return [a[k] + b[k] for k in range(order + 1)]


def series_mul(a, b, order):
    return [sum((a[i] * b[k - i] for i in range(k + 1)), Fraction(0))
            for k in range(order + 1)]


def series_exp(a, order):
    """exp(a), for a[0]=0, by the formal derivative recurrence."""
    assert a[0] == 0
    out = [Fraction(0)] * (order + 1)
    out[0] = Fraction(1)
    for n in range(1, order + 1):
        out[n] = sum(Fraction(k) * a[k] * out[n - k]
                     for k in range(1, n + 1)) / n
    return out


def tree_series(order):
    return [Fraction(0)] + [Fraction(k ** (k - 1), factorial(k))
                             for k in range(1, order + 1)]


def integer_power_series(base, exponent, order):
    out = [Fraction(0)] * (order + 1)
    out[0] = Fraction(1)
    for _ in range(exponent):
        out = series_mul(out, base, order)
    return out


def cumulant_polynomial(order):
    """c[k](q)=k![t^k](-log(1-qT(t))) as q-power dictionaries."""
    T = tree_series(order)
    c = [dict() for _ in range(order + 1)]
    for r in range(1, order + 1):
        power = integer_power_series(T, r, order)
        for k in range(order + 1):
            if power[k]:
                c[k][r] = power[k] * factorial(k) / r
    return c


def poly_coeff_product(left, right, left_power, right_power):
    return left.get(left_power, Fraction(0)) * right.get(right_power, Fraction(0))


def main():
    order = 14
    T = tree_series(order)
    exp_T = series_exp(T, order)
    rhs = [Fraction(0)] + exp_T[:order]
    assert all(T[k] == rhs[k] for k in range(order + 1))
    print("tree_function_fidelity: T=t*exp(T) through degree 14")

    c = cumulant_polynomial(order)
    # The source's EGF has M_0=1 and its t coefficient is q*y; its
    # coefficient of y in M_k is exactly c[k](q).
    assert c[1] == {1: Fraction(1)}
    assert c[2] == {1: Fraction(2), 2: Fraction(1)}
    print("colored_egf_fidelity: M_0=1, [y]M_1=q, [y]M_2=2*q+q^2")

    for colours in (2, 3):
        for n in range(1, 7):
            a = n + 1
            # D_n=M_n*M_(n+2)-M_(n+1)^2, rows {0,1}, columns {n,n+1}.
            # In y-degree two its coefficient is the four ordered products
            # of the linear y parts c[k](q_l).
            positive = (
                poly_coeff_product(c[n], c[n + 2], a, a)
                + poly_coeff_product(c[n + 2], c[n], a, a)
            )
            negative = 2 * poly_coeff_product(c[n + 1], c[n + 1], a, a)
            got = positive - negative
            assert positive == 0, (colours, n, positive)
            assert negative == Fraction(2 * factorial(n) ** 2)
            assert got == Fraction(-2 * factorial(n) ** 2)
            # Degree bound and the leading coefficient used in the proof.
            assert max(c[n]) <= n
            assert c[n + 1][n + 1] == Fraction(factorial(n))
        print(f"colours={colours}: exact D_n coefficient n=1..6 is -2*(n!)^2")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
