#!/usr/bin/env python3
"""Exact rational checks for Moroz arXiv:2510.03373v2, Sections 2 and 4.1.

The source alternating summand ends in (q_k - 1). All arithmetic below
uses integers and Fraction. Finite tests support, but do not prove, the
unbounded differentiability assertion.
"""

from fractions import Fraction as Q
from itertools import product
from math import prod


def cylinder(digits, parameters):
    """Positive partial sum, alternating partial sum, and common diameter."""
    assert len(digits) == len(parameters)
    plus, minus, length = Q(0), Q(0), Q(1)
    for k, (d, r) in enumerate(zip(digits, parameters)):
        assert 1 <= r < d
        plus += length * Q(r, d)
        minus += (-1) ** k * length * Q(r, d - 1)
        length *= Q(r, (d - 1) * d)
    plus_interval = (plus, plus + length)
    end = minus + (-1) ** len(digits) * length
    minus_interval = (min(minus, end), max(minus, end))
    return plus_interval, minus_interval, length, minus


def source_terms(digits, parameters):
    """Direct products, independent of the cylinder recurrence."""
    plus, minus = Q(0), Q(0)
    for k, d in enumerate(digits):
        numer = prod(parameters[: k + 1])
        denom = prod((a - 1) * a for a in digits[:k])
        plus += Q(numer, denom * d)
        minus += (-1) ** k * Q(numer, denom * (d - 1))
    return plus, minus


def periodic_values(d):
    """P = 1, constant digit d; exact sums of geometric series."""
    ratio = Q(1, d * (d - 1))
    return Q(1, d) / (1 - ratio), Q(1, d - 1) / (1 + ratio)


def eventually_constant(prefix, tail):
    ip, _, length, a = cylinder(prefix, [1] * len(prefix))
    x, fx = periodic_values(tail)
    return ip[0] + length * x, a + (-1) ** len(prefix) * length * fx


def fidelity():
    assert periodic_values(2) == (Q(1), Q(2, 3))
    assert periodic_values(3) == (Q(2, 5), Q(3, 7))
    assert periodic_values(5) == (Q(4, 19), Q(5, 21))
    count = 0
    for n in range(1, 5):
        for rows in product(product(range(1, 4), range(3)), repeat=n):
            rs = [r for r, _ in rows]
            ds = [r + 1 + off for r, off in rows]
            ip, im, length, a = cylinder(ds, rs)
            direct_plus, direct_minus = source_terms(ds, rs)
            diameter = Q(prod(rs), prod((d - 1) * d for d in ds))
            assert ip[0] == direct_plus and a == direct_minus
            assert ip[1] - ip[0] == im[1] - im[0] == diameter == length
            assert 0 < length <= Q(1, 2**n)
            count += 1
    # The displayed ratio in the source's Section 3 is reciprocated.
    l2 = cylinder([2], [1])[2]
    l3 = cylinder([3], [1])[2]
    assert l2 / l3 == 3 and l2 / l3 != Q(1, 3)
    print(f"FIDELITY: {count} source-product and equal-diameter cases; shrink <= 2^-n")
    print("Luroth: constant digits 2,3,5 give (1,2/3), (2/5,3/7), (4/19,5/21)")
    print("SOURCE TYPO CONTROL: |Delta_2|/|Delta_3| = 3, not 1/3")


def local_secants():
    count = 0
    worst_ratio = Q(0)
    minimum_signed_secant = None
    for n in range(1, 5):
        for rows in product(product(range(1, 4), range(3)), repeat=n):
            rs = [r for r, _ in rows]
            ds = [r + 1 + off for r, off in rows]
            ip, im, li, _ = cylinder(ds, rs)
            js = ds[:-1] + [ds[-1] + 1]
            jp, jm, lj, a = cylinder(js, rs)
            ratio = li / lj
            assert ratio == Q(ds[-1] + 1, ds[-1] - 1) <= 3
            worst_ratio = max(worst_ratio, ratio)
            assert jp[1] == ip[0]
            assert (jm[1] == im[0]) if n % 2 else (im[1] == jm[0])
            # Complete the code of J with r = 1 and digit 2 forever.
            # Its exact positive value is J's right endpoint; its alternating
            # value is 2/3 of J's diameter toward the required far endpoint.
            y, fy = jp[1], a + (-1) ** n * lj * Q(2, 3)
            assert jm[0] < fy < jm[1]
            # The estimate uses only membership, so independent coordinates
            # here cover every pair satisfying these geometric constraints.
            for tx, tf in product([Q(1, 7), Q(6, 7), Q(1)], [Q(1, 7), Q(6, 7)]):
                x = ip[0] + tx * li
                fx = im[0] + tf * li
                assert 0 < x - y <= li + lj <= 4 * lj
                quotient = (fx - fy) / (x - y)
                signed = quotient if n % 2 else -quotient
                assert signed > Q(1, 8)
                minimum_signed_secant = signed if minimum_signed_secant is None else min(
                    minimum_signed_secant, signed
                )
                count += 1
    print(f"LOCAL SECANTS: {count} exact geometric cases; largest adjacent ratio = {worst_ratio}")
    print(f"LOCAL SECANTS: minimum signed sampled quotient = {minimum_signed_secant} > 1/8")


def luroth_secants():
    cases = [([], 2), ([], 3), ([], 5), ([2], 3), ([3], 2)]
    count = 0
    for prefix, tail in cases:
        x, fx = eventually_constant(prefix, tail)
        residual, digits, parameters = x, [], []
        first = []
        for n in range(1, 17):
            inv = 1 / residual
            d = inv.numerator // inv.denominator + 1
            digits.append(d)
            parameters.append(1)
            residual = (residual - Q(1, d)) * d * (d - 1)
            assert 0 < residual <= 1
            ip, im, li, _ = cylinder(digits, parameters)
            assert ip[0] < x <= ip[1] and im[0] < fx < im[1]
            js = digits[:-1] + [d + 1]
            jp, _, lj, a = cylinder(js, parameters)
            y = jp[1]
            fy = a + (-1) ** n * lj * Q(2, 3)
            # Re-evaluate the actual eventually-constant code independently.
            assert (y, fy) == eventually_constant(js, 2)
            assert 0 < x - y <= Q(1, 2 ** (n - 1))
            quotient = (fx - fy) / (x - y)
            assert quotient > Q(1, 8) if n % 2 else quotient < -Q(1, 8)
            if n <= 4:
                first.append(str(quotient))
            count += 1
        print(f"Luroth x={x}, F(x)={fx}: first four left secants = {', '.join(first)}")
    print(f"LUROTH SECANTS: {count} actual exact pairs, depths 1..16")


if __name__ == "__main__":
    fidelity()
    local_secants()
    luroth_secants()
    print("ALL_CHECKS_OK")
