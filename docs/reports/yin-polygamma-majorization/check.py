#!/usr/bin/env python3
"""Exact checks for Yin, JMI 16(3) (2022), Conjecture 4.1.

The derivative definition and equation (1.5) bind the functions. The conflicting
finite-sum sign in equation (1.6) is diagnosed, not used as a premise.
"""
from fractions import Fraction as F
from math import factorial, prod
import sympy as sp


def moment(p, k, x, r):
    return sum((F(1) / (x + j * k) ** r for j in range(p + 1)), F(0))


def derivative(p, k, x, r):
    assert p >= 1 and k > 0 and x > 0 and r >= 1
    return (-1) ** (r + 1) * factorial(r) * moment(p, k, x, r + 1)


def quotient(p, k, x, alpha, beta):
    return prod(derivative(p, k, x, r) for r in alpha) / prod(
        derivative(p, k, x, r) for r in beta)


def strictly_majorized(alpha, beta):
    a, b = sorted(alpha, reverse=True), sorted(beta, reverse=True)
    return (a != b and sum(a) == sum(b)
            and all(sum(a[:r]) <= sum(b[:r]) for r in range(1, len(a))))


def main():
    x, k = sp.symbols('x k', positive=True)
    fidelity_count = 0
    bridge_count = 0
    for p in range(1, 5):
        # Source (1.4), including its normalization and exponent x/k - 1.
        gamma = (factorial(p + 1) * k ** (p + 1) * (p * k) ** (x / k - 1)
                 / sp.prod(x + j * k for j in range(p + 1)))
        psi = sp.log(p * k) / k - sum(1 / (x + j * k) for j in range(p + 1))
        assert sp.simplify(sp.diff(gamma, x) / gamma - psi) == 0
        fidelity_count += 1
        for r in range(1, 7):
            rhs = (-1) ** (r + 1) * factorial(r) * sum(
                1 / (x + j * k) ** (r + 1) for j in range(p + 1))
            assert sp.simplify(sp.diff(psi, x, r) - rhs) == 0
            bridge_count += 1
    print(f'FIDELITY_OK: source (1.4) logarithmic derivative equals (1.5), {fidelity_count} symbolic p cases')
    print(f'DERIVATIVE_BRIDGE_OK: {bridge_count} symbolic (p,r) cases')
    assert derivative(1, F(1), F(1), 1) == F(5, 4)
    assert -factorial(1) * moment(1, F(1), F(1), 2) == F(-5, 4)
    print('PRINTED_SIGN_CONFLICT_CONFIRMED: (1.6) finite sum gives -5/4; derivative definition gives 5/4 at p=k=x=r=1')

    a, b = (2, 2, 2), (1, 2, 3)
    assert strictly_majorized(a, b)
    q1, q2 = quotient(1, F(1), F(1), a, b), quotient(1, F(1), F(2), a, b)
    assert q1 == F(54, 85) and q2 == F(2450, 3783)
    assert q2 - q1 == F(3968, 321555) > 0
    print(f'COUNTEREXAMPLE_OK: Q(1)={q1}; Q(2)={q2}; Q(2)-Q(1)={q2-q1}')

    count = 0
    for p in range(1, 13):
        for m in range(2, 13):
            for kv in (F(1, 3), F(1), F(3, 2)):
                for xv in (F(1, 7), F(1), F(2), F(11)):
                    u = [F(1) / (xv + j * kv) for j in range(p + 1)]
                    sm, sn, so = (sum(v ** r for v in u) for r in (m, m + 1, m + 2))
                    gap = sm * so - sn ** 2
                    pair_gap = sum(u[i] ** m * u[j] ** m * (u[i] - u[j]) ** 2
                                   for i in range(p + 1) for j in range(i + 1, p + 1))
                    assert gap == pair_gap > 0
                    ratio = quotient(p, kv, xv, (m, m), (m - 1, m + 1))
                    assert ratio == F(m, m + 1) * sn ** 2 / (sm * so)
                    assert 0 < ratio < F(m, m + 1)
                    count += 1
    assert count == 1584
    print(f'MOMENT_GAP_AND_STRICT_DEFICIT_OK: {count} exact rational (p,m,k,x) cases')

    embed_count = 0
    for n in range(3, 9):
        for m in range(2, 9):
            a = (m,) * n
            b = (m - 1, m + 1) + (m,) * (n - 2)
            assert strictly_majorized(a, b)
            for p in range(1, 5):
                assert quotient(p, F(1), F(2), a, b) == quotient(
                    p, F(1), F(2), (m, m), (m - 1, m + 1))
                embed_count += 1
    print(f'MULTIFACTOR_EMBEDDING_OK: {embed_count} exact cases, n=3..8')

    endpoint_count = 0
    for p in range(1, 4):
        for m in range(2, 5):
            s = lambda r: sum(1 / (x + j) ** r for j in range(p + 1))
            ratio = sp.cancel(sp.Rational(m, m + 1) * s(m + 1) ** 2 / (s(m) * s(m + 2)))
            assert sp.limit(ratio, x, 0, dir='+') == sp.Rational(m, m + 1)
            assert sp.limit(ratio, x, sp.oo) == sp.Rational(m, m + 1)
            endpoint_count += 1
    print(f'ENDPOINTS_OK: {endpoint_count} exact rational-function cases, both limits m/(m+1)')
    print('SCOPE: finite checks are regression evidence; the unbounded family requires the general argument in README and preregistration')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
