#!/usr/bin/env python3
"""Exact checks for the q-Narayana strong-log-convexity question."""

from fractions import Fraction
from math import factorial

import sympy as sp


def qpoch(a: Fraction, n: int) -> Fraction:
    out = Fraction(1)
    for j in range(n):
        out *= 1 - a * (q := Fraction(1, 2)) ** j
    return out


def qpoch_at(a: Fraction, q: Fraction, n: int) -> Fraction:
    out = Fraction(1)
    for j in range(n):
        out *= 1 - a * q**j
    return out


def qbinom(n: int, k: int, q: Fraction) -> Fraction:
    if k < 0 or k > n:
        return Fraction(0)
    return qpoch_at(q, q, n) / (qpoch_at(q, q, k) * qpoch_at(q, q, n - k))


def kernel(n: int, k: int, theta: int, q: Fraction) -> Fraction:
    return (
        q ** (k * (k - 1))
        * qbinom(n, k, q)
        * qpoch_at(q ** (theta + 1), q, n)
        / (qpoch_at(q ** (theta + 1), q, k) * qpoch_at(q ** (theta + 1), q, n - k))
    )


def rows(theta: int, q: Fraction, upto: int):
    return [[kernel(n, k, theta, q) for k in range(n + 1)] for n in range(upto + 1)]


def transform(row, inputs):
    return sum((x * a for x, a in zip(row, inputs)), Fraction(0))


def main() -> None:
    q = Fraction(1, 2)
    r = rows(theta=1, q=q, upto=3)
    assert r == [
        [Fraction(1)],
        [Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(7, 4), Fraction(1, 4)],
        [Fraction(1), Fraction(35, 16), Fraction(35, 64), Fraction(1, 64)],
    ], r

    inputs = [factorial(n + 1) for n in range(4)]
    outputs = [transform(row, inputs[: len(row)]) for row in r]
    assert outputs == [Fraction(1), Fraction(3), Fraction(6), Fraction(289, 32)]
    assert outputs[0] * outputs[2] - outputs[1] ** 2 == -3
    assert outputs[1] * outputs[3] - outputs[2] ** 2 == Fraction(-285, 32)

    qsym = sp.symbols("q")
    b2 = 1 + (1 + qsym + qsym**2) * 2 + qsym**2 * 6
    delta = sp.expand(b2 - 3**2)
    assert delta == 8 * qsym**2 + 2 * qsym - 6
    assert sp.Poly(delta, qsym).coeff_monomial(1) == -6
    assert sp.expand(sp.factor(delta) - 2 * (qsym + 1) * (4 * qsym - 3)) == 0

    # Structural obstruction: rows [1], [1,1], [1,C,d] and geometric input.
    for theta in (0, 1, 2, 3):
        for q0 in (Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
            c = (1 + q0) * (1 - q0 ** (theta + 2)) / (1 - q0 ** (theta + 1))
            d = q0**2
            assert c > 0 and d < 1
            scale = (c + 1) / (1 - d)
            b0 = Fraction(1)
            b1 = 1 + scale
            b2s = 1 + c * scale + d * scale**2
            assert b0 * b2s - b1**2 == -3 * scale

    # Strict strong log-convexity of the factorial input.
    for i in range(1, 5):
        for j in range(i, 5):
            lhs = factorial(i) * factorial(j + 2)
            rhs = factorial(i + 1) * factorial(j + 1)
            assert lhs - rhs == factorial(i) * factorial(j + 1) * (j - i + 1) > 0

    print("theta=1, q=1/2 rows:", r)
    print("factorial outputs:", outputs)
    print("formal q minor:", delta)
    print("structural family and strict input checks: ok")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
