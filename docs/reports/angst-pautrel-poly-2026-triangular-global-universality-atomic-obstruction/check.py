#!/usr/bin/env python3
"""Exact checks for Angst--Pautrel--Poly, arXiv:2609.01007v1, Remark 5(2).

Input source: APP_revision.tex, §§1.1--1.3. All arithmetic used as evidence
is integer, Fraction, or exact SymPy algebra. No simulation or floating point.
Finite checks support the proposed outline; they do not prove its unbounded claim.
"""
from fractions import Fraction as Q
from itertools import product
import sympy as sp


def triangle_pi_units(q):
    """f(pi*q), reducing q to [-1,1); endpoints both have value -1."""
    q = Q(q)
    r = (q + 1) % 2 - 1
    return 1 + 2 * r if r <= 0 else 1 - 2 * r


def signal_pi_units(coefficients, q):
    return sum((a * triangle_pi_units(k * q)
                for k, a in enumerate(coefficients, 1)), Q(0))


def check_fidelity():
    # Source §1.2: the continuous periodic triangle interpolates cos on pi*Z.
    for k in range(-128, 129):
        assert triangle_pi_units(k) == (1 if k % 2 == 0 else -1)
    for denominator in range(1, 18):
        for numerator in range(-3 * denominator, 3 * denominator + 1):
            q = Q(numerator, denominator)
            assert triangle_pi_units(q + 2) == triangle_pi_units(q)
            if -1 <= q <= 0:
                assert triangle_pi_units(q) == 1 + 2 * q
            if 0 <= q <= 1:
                assert triangle_pi_units(q) == 1 - 2 * q
    # Source Remark 5(1): cosine gives the constant 2/sqrt(3).
    x = sp.Symbol('x', real=True)
    cosine_norm = sp.integrate(sp.cos(x) ** 2, (x, 0, 2 * sp.pi))
    derivative_norm = sp.integrate(sp.diff(sp.cos(x), x) ** 2,
                                   (x, 0, 2 * sp.pi))
    assert cosine_norm == derivative_norm == sp.pi
    cosine_constant = 2 / sp.sqrt(3) * sp.sqrt(derivative_norm / cosine_norm)
    assert sp.simplify(cosine_constant - 2 / sp.sqrt(3)) == 0
    print('FIDELITY_OK: periodic triangle interpolates cosine on pi*Z; '
          'Remark 5(1) cosine constant = 2/sqrt(3)')

    left = 1 + 2 * x / sp.pi
    right = 1 - 2 * x / sp.pi
    mean = sp.integrate(left, (x, -sp.pi, 0)) + sp.integrate(right, (x, 0, sp.pi))
    l2 = sp.integrate(left ** 2, (x, -sp.pi, 0)) + sp.integrate(right ** 2, (x, 0, sp.pi))
    d2 = (sp.integrate(sp.diff(left, x) ** 2, (x, -sp.pi, 0))
          + sp.integrate(sp.diff(right, x) ** 2, (x, 0, sp.pi)))
    assert mean == 0 and l2 == 2 * sp.pi / 3 and d2 == 8 / sp.pi
    assert sp.simplify(2 / sp.sqrt(3) * sp.sqrt(d2 / l2) - 4 / sp.pi) == 0
    print('TRIANGLE_NORMALIZATION_OK: integral f = 0; <f,f> = 2*pi/3; '
          "<f',f'> = 8/pi; proposed constant = 4/pi")


def check_family():
    j, u, c = sp.symbols('j u c', real=True)
    indices = [4 * j + r for r in range(1, 5)]
    signs = [1, -1, -1, 1]
    assert sum(signs) == 0
    assert sp.expand(sum(a * k for a, k in zip(signs, indices))) == 0
    assert sp.expand(sum(c * a * (1 - 2 * k * u)
                         for a, k in zip(signs, indices))) == 0
    assert [triangle_pi_units(Q(k, 2)) for k in range(1, 5)] == [0, -1, 0, 1]
    print('SYMBOLIC_BLOCK_OK: moments (0,0); affine block identically 0; '
          'quarter-period block value = 2*c')

    for m in range(1, 101):
        n = 4 * m
        for amplitude in [Q(1), Q(2, 3), Q(7, 5)]:
            coefficients = [amplitude * a for a in signs] * m
            assert all(a != 0 for a in coefficients)
            assert sum(coefficients) == 0
            assert sum(k * a for k, a in enumerate(coefficients, 1)) == 0
            for q in [Q(0), Q(1, 4 * n), Q(1, 2 * n), Q(3, 4 * n), Q(1, n)]:
                assert 0 <= q <= Q(1, n)
                assert signal_pi_units(coefficients, q) == 0
            assert signal_pi_units(coefficients, Q(1, 2)) == 2 * amplitude * m
        prescribed_probability = Q(1, 2) ** n
        assert prescribed_probability > 0
        for epsilon in [Q(1, 2), Q(1, 10), Q(1, 1000)]:
            assert (epsilon / 2) ** n > 0
        if m in [1, 2, 3, 10, 100]:
            print(f'FAMILY_CASE_OK: m={m}, n={n}, S_n(pi/2)={2*m}, '
                  f'prescribed Rademacher event probability=1/2^{n}')
    assert (Q(1, 2) * 1 + Q(1, 2) * -1) == 0
    assert (Q(1, 2) * 1 ** 2 + Q(1, 2) * (-1) ** 2) == 1
    assert (Q(1, 2) * abs(1) ** 3 + Q(1, 2) * abs(-1) ** 3) == 1
    print('FAMILY_RANGE_OK: m=1..100, 3 nonzero amplitudes, 5 exact interval '
          'samples each; Rademacher moments (mean,second,third absolute)=(0,1,1)')

    # Independent enumeration of all signs, without assuming the prescribed pattern.
    for n in [4, 8, 12]:
        flat = 0
        for coefficients in product([-1, 1], repeat=n):
            if (sum(coefficients) == 0
                    and sum(k * a for k, a in enumerate(coefficients, 1)) == 0):
                flat += 1
        assert flat >= 2
        print(f'ENUMERATION_OK: n={n}, flat-first-cell sign vectors={flat}, '
              f'probability={Q(flat, 2**n)}')


if __name__ == '__main__':
    check_fidelity()
    check_family()
    print('ALL_CHECKS_OK')
