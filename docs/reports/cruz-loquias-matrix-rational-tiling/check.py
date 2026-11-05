#!/usr/bin/env python3
"""Exact finite checks for the fixed-root Cruz--Loquias v1 tiling candidate.

This is experimental evidence for the preregistered argument, not a Lean proof.
Only integers and fractions.Fraction are used.
"""
from fractions import Fraction as F
from math import gcd


def matvec(a, x):
    return tuple(sum((F(c) * y for c, y in zip(row, x)), F(0)) for row in a)


def fidelity():
    # v1 Section 4.2 example: the path label 443 ends at (-5,3).
    p = ((4, -1), (1, 1))
    qinv = ((F(1, 2), F(-5, 2)), (0, 1))
    x = (F(0), F(0))
    vertices = [x]
    for a in (4, 4, 3):
        y = matvec(p, x)
        x = matvec(qinv, (y[0] + a, y[1]))
        assert all(t.denominator == 1 for t in x)
        vertices.append(x)
    assert x == (-5, 3)
    print('FIDELITY v1 Section 4.2 path 443:', vertices)


def phi(x, m):
    y = 2 * x
    a = (y.numerator * pow(y.denominator, -1, m)) % m
    return (y - a) / -m


def children(n, m):
    return {(-m * n + a) // 2 for a in range(m) if (-m * n + a) % 2 == 0}


def ceilq(x):
    return -((-x.numerator) // x.denominator)


def family(m):
    assert m >= 9 and m % 2 == 1 and gcd(m, 2) == 1
    b, radius = F(m, 2), F(m - 1, m - 2)
    # Bezout: (-m)*1 + 2*((m+1)/2) = 1, in every dimension.
    assert -m + 2 * ((m + 1) // 2) == 1 and b > 1
    for r in range(7):
        for t in range(-50, 51):
            x = F(t, 2**r)
            if x.denominator > 1:
                assert phi(x, m).denominator <= x.denominator // 2
            seen = set()
            while x:
                assert x not in seen
                seen.add(x)
                x = phi(x, m)
            assert len(seen) <= r + 12
    for n in range(-300, 301):
        assert phi(F(n), m) == -(2 * n // m)
        if n >= 1 or n <= -2:
            assert abs(phi(F(n), m)) < abs(n)
    assert phi(F(-1), m) == 1 and phi(F(1), m) == 0
    l = u = 0
    states = {0}
    certificates = []
    for j in range(1, 13):
        old_l, old_u = l, u
        l = ceilq(-b * old_u)
        u = (-b * old_l + F(m - 1, 2)).__floor__()
        assert u - l >= b * (old_u - old_l) + F(m - 5, 2)
        if j <= 3:
            states = set().union(*(children(n, m) for n in states))
            assert states == set(range(l, u + 1))
        if j == 1:
            assert u - l == (m - 1) // 2
        # If the preregistered interval/limit lemma holds, this is an exact
        # lower bound on its width, from the finite hull and uniform tail error.
        lower = F(u - l, 1) / b**j - 2 * radius / b**j
        if lower > 1:
            certificates.append((j, lower))
    assert certificates
    bound = 1 + F(m - 8, m * (m - 2))
    summed = F(m - 1, m) + F(2 * (m - 5), m * (m - 2))
    assert bound == summed and bound > 1
    for d in range(2, 6):
        assert (bound - 1) * bound**(d - 1) > 0
        assert (b * bound - 1) * (b * bound)**(d - 1) > 0
    # Distinct one-digit integer translations, and admissibility from the root.
    assert 2 in range(m) and 4 in range(m)
    assert F(2, 2) == 1 and F(4, 2) == 2
    assert 1 in children(0, m) and 2 in children(0, m)
    print(f'FAMILY m={m}: width_bound={bound}; first_finite_tail_certificate={certificates[0]}; d=2..5 overlap_bounds>0; printed_power_zero_scale={b}')


def parity_identity_checks():
    # Floor/ceiling error depends only on endpoint parity. Testing the four
    # cases plus the displayed exact affine formula checks the recurrence.
    for m in range(9, 102, 2):
        c = F(m - 1, 2)
        for el in (0, 1):
            for eu in (0, 1):
                l, u = -20 + el, 20 + eu
                actual = (-F(m, 2) * l + c).__floor__() - ceilq(-F(m, 2) * u)
                error = actual - F(m, 2) * (u - l)
                formula = c - (F(1, 2) if el else 0) - (F(1, 2) if eu else 0)
                assert error == formula and error >= F(m - 5, 2)
    print('PARITY recurrence: all four endpoint parities, odd m=9..101, exact floor-error formula verified')


if __name__ == '__main__':
    fidelity()
    parity_identity_checks()
    for m in range(9, 50, 2):
        family(m)
    print('SCOPE finite checks only: odd m=9..49, depths<=12, dyadic r<=6 |t|<=50, integers |n|<=300; no floating point')
    print('ALL_CHECKS_OK')
