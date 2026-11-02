#!/usr/bin/env python3
"""Exact checks for Li's k-factor Kronecker separation candidate.

All quantities are integer exponents of 2.  No floating point is used.
"""
from fractions import Fraction
from itertools import product
import sys


def h3(z, u, v):
    z1, z2, z3 = z
    return 3 * u * (z1 * z2 + z1 * z3 + z2 * z3) + 3 * v * z1 * z2 * z3


def all_exponents(k, u, v):
    # H is the 3-factor witness tensored with identity factors.
    for z in product((-1, 1), repeat=k):
        yield z, h3(z[:3], u, v)


def width(values):
    return max(values) - min(values)


def pow2(exp):
    """Return 2**exp exactly, including for negative exponents."""
    return 2**exp if exp >= 0 else Fraction(1, 2**(-exp))


def check_explicit_three_factor():
    # The source's concrete example is H = diag(2^12,2^-6,2^-6,1,2^-6,1,1,2^6).
    h = [12, -6, -6, 0, -6, 0, 0, 6]
    assert h == [h3(z, 1, 1) for z in product((1, -1), repeat=3)]
    assert width(h) == 18
    # G* = diag(2^3,2,2,2^-1,2,2^-1,2^-1,2^-3).
    g = [3, 1, 1, -1, 1, -1, -1, -3]
    relative = [x - y for x, y in zip(h, g)]
    assert relative == [9, -7, -7, 1, -7, 1, 1, 9]
    assert width(relative) == 16
    assert 2**width(h) == 262144
    assert 2**width(relative) == 65536
    assert [pow2(x) for x in relative] == [512, Fraction(1, 128), Fraction(1, 128), 2,
                                         Fraction(1, 128), 2, 2, 512]


def check_nearest_identity(u, v, k):
    # Every one-mode partial trace of log_2(H) vanishes.
    values = list(all_exponents(k, u, v))
    for j in range(k):
        for sign in (-1, 1):
            partial = sum(value for z, value in values if z[j] == sign)
            assert partial == 0, (k, j, sign, partial)


def check_exact_additive_optimum(u, v, k):
    # The residual is h - sum_j a_j(z_j).  The exact dual witness is the
    # three-factor signed measure, averaged uniformly over extra coordinates.
    values = list(all_exponents(k, u, v))
    witness = {}
    base_w = {
        (1, 1, 1): Fraction(2, 3),
        (-1, -1, -1): Fraction(1, 3),
        (-1, 1, 1): Fraction(-1, 3),
        (1, -1, 1): Fraction(-1, 3),
        (1, 1, -1): Fraction(-1, 3),
    }
    for z, _ in values:
        witness[z] = base_w.get(z[:3], Fraction(0)) / (2 ** (k - 3))
    assert sum(witness.values()) == 0
    assert sum(x for x in witness.values() if x > 0) == 1
    assert sum(-x for x in witness.values() if x < 0) == 1
    for j in range(k):
        for sign in (-1, 1):
            assert sum(w for z, w in witness.items() if z[j] == sign) == 0

    # Additive factors cancel under the witness.  Therefore every residual
    # range is at least witness dot h = 12u+4v.
    dot_h = sum(witness[z] * h for z, h in values)
    assert dot_h == 12 * u + 4 * v

    # The proposed local correction is a_j(z_j) = v z_j.  Its residual width
    # is computed exactly and reaches the lower bound when v < 2u.
    residual = [h - v * sum(z[:3]) for z, h in values]
    assert width(residual) == 12 * u + 4 * v, (k, u, v, width(residual))
    assert min(residual) == -3 * u - 4 * v
    assert max(residual) == 9 * u


def main():
    check_explicit_three_factor()
    for u, v in ((1, 1), (2, 1), (3, 1), (3, 5), (7, 10)):
        assert v < 2 * u
        for k in range(3, 7):
            check_nearest_identity(u, v, k)
            check_exact_additive_optimum(u, v, k)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"CHECK_FAILED: {exc}", file=sys.stderr)
        raise
