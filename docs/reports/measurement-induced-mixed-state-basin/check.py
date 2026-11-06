#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact fidelity and Lyapunov checks for arXiv:1808.08151v1, Eq. (10).

No floating-point values or numerical convergence tests are evidence here.
Symbolic rational-function identities are supplemented by finite rational tests;
neither is a Lean proof of the universally quantified convergence claim.
"""
from fractions import Fraction as Q
from itertools import product
import sympy as s


def zero(expr):
    result = s.cancel(s.together(expr))
    assert result == 0, result


def matrix_zero(matrix):
    for entry in matrix:
        zero(entry)


def nonnegative_coefficients(expr, variables):
    coefficients = s.Poly(s.expand(expr), *variables).coeffs()
    assert all(coefficient >= 0 for coefficient in coefficients)


def step(r):
    u, v, w = r
    h = 1 + w * w
    assert h > 0
    return ((u * u - v * v) / h, 2 * w / h, -2 * u * v / h)


def norm_sq(r):
    return sum((t * t for t in r), Q(0))


def transformed(r):
    d = 1 - norm_sq(r)
    assert d > 0
    return tuple(t * t / d for t in r)


def tstep(t):
    A, B, C = t
    D = 1 + 2 * A + 2 * B
    return ((A - B) ** 2 / D, 4 * C * (1 + A + B + C) / D,
            4 * A * B / D)


def V(t):
    A, B, C = t
    return A * B + A * C + B * C + C * C + A / 2 + B / 4 + C


def residual(t):
    A, B, C = t
    return ((A + B) ** 2 * (A + 2 * C + Q(1, 2))
            + 2 * C * (A + B) * (C + 1) + A * A + A / 2 + B / 4)


def symbolic_checks():
    u, v, w = s.symbols('u v w', real=True)
    rho = s.Matrix([[1 + w, u - s.I * v], [u + s.I * v, 1 - w]]) / 2
    h = s.expand(rho[0, 0] ** 2 + rho[1, 1] ** 2)
    zero(h - (1 + w ** 2) / 2)
    squared = rho.applyfunc(lambda entry: entry ** 2) / h
    U = s.Matrix([[1, s.I], [s.I, 1]]) / s.sqrt(2)
    matrix_zero(U * U.conjugate().T - s.eye(2))
    final = U * squared * U.conjugate().T
    # Read output coordinates from the density matrix; Eq. (10) is the target.
    extracted = (s.expand(final[0, 1] + final[1, 0]),
                 s.expand(s.I * (final[0, 1] - final[1, 0])),
                 s.expand(final[0, 0] - final[1, 1]))
    expected = ((u ** 2 - v ** 2) / (1 + w ** 2),
                2 * w / (1 + w ** 2), -2 * u * v / (1 + w ** 2))
    for actual, target in zip(extracted, expected):
        zero(actual - target)
    target_rho = s.Matrix([[1 + expected[2], expected[0] - s.I * expected[1]],
                           [expected[0] + s.I * expected[1], 1 - expected[2]]]) / 2
    matrix_zero(final - target_rho)
    zero(s.trace(rho * rho) - (1 + u ** 2 + v ** 2 + w ** 2) / 2)
    zero(rho.det() - (1 - u ** 2 - v ** 2 - w ** 2) / 4)
    print('FIDELITY: protocol derives all three coordinates of source Eq. (10); residuals 0')
    print('FIDELITY: source Eq. (9) purity and determinant; residuals 0')

    a, b, c = s.symbols('a b c')
    d = 1 - a - b - c
    ap, bp, cp = (a - b) ** 2 / (1 + c) ** 2, 4 * c / (1 + c) ** 2, 4 * a * b / (1 + c) ** 2
    dp = 1 - ap - bp - cp
    zero(dp - d * (d + 2 * a + 2 * b) / (1 + c) ** 2)
    print('DEFECT: d_next = d (d + 2a + 2b)/(1+c)^2; residual 0')
    A, B, C = s.symbols('A B C')
    inverse = {a: A / (1 + A + B + C), b: B / (1 + A + B + C),
               c: C / (1 + A + B + C)}
    tp = tstep((A, B, C))
    for q, target in zip((ap, bp, cp), tp):
        zero(s.cancel(q / dp).subs(inverse, simultaneous=True) - target)
    print('TRANSFORM: all three rational coordinate conjugacies; residuals 0')
    D = 1 + 2 * A + 2 * B
    R = residual((A, B, C))
    zero(D ** 2 * (V((A, B, C)) - V(tp)) - R)
    print('LYAPUNOV: D^2 (V - V_next) = displayed nonnegative certificate; residual 0')
    for expr in (R, R - (A + B) / 4, V((A, B, C)) - (A + B) / 4,
                 V((A, B, C)) - C, 8 * V((A, B, C)) - 2 * (A + B)):
        nonnegative_coefficients(expr, (A, B, C))
    zero(D * tp[1] - 4 * C - 4 * C * (A + B + C))
    nonnegative_coefficients(4 * C * (A + B + C), (A, B, C))
    print('COERCIVITY: R >= (A+B)/4; V >= (A+B)/4, C; D <= 1+8V; exact certificates')
    print('C-CONTROL: D B_next - 4C = 4C(A+B+C) >= 0; exact certificate')


def rational_checks():
    values = set()
    for denominator in range(1, 7):
        for nums in product(range(-denominator, denominator + 1), repeat=3):
            r = tuple(Q(k, denominator) for k in nums)
            if norm_sq(r) < 1:
                values.add(r)
    for r in sorted(values):
        rp = step(r)
        d, dp = 1 - norm_sq(r), 1 - norm_sq(rp)
        a, b, c = (t * t for t in r)
        assert dp == d * (d + 2 * a + 2 * b) / (1 + c) ** 2 > 0
        assert transformed(rp) == tstep(transformed(r))
    print(f'RATIONAL FIDELITY: {len(values)} distinct interior triples, common denominators 1..6')

    grid = (Q(0), Q(1, 4), Q(1, 2), Q(1), Q(2), Q(4))
    orbit_count = 0
    for initial in product(grid, repeat=3):
        t = initial
        M = V(initial)
        K = 1 + 8 * M
        partial = Q(0)
        for _ in range(6):
            A, B, C = t
            D = 1 + 2 * A + 2 * B
            tn = tstep(t)
            assert all(z >= 0 for z in t)
            assert 0 <= V(tn) <= V(t) <= M
            assert A + B <= 4 * M and C <= M and 1 <= D <= K
            assert D ** 2 * (V(t) - V(tn)) == residual(t)
            assert V(t) - V(tn) >= (A + B) / (4 * K ** 2)
            assert C <= K * tn[1] / 4
            partial += A + B
            assert partial <= 4 * K ** 2 * (M - V(tn))
            t = tn
        orbit_count += 1
    print(f'RATIONAL LYAPUNOV: {orbit_count} nonnegative initial triples, six exact steps each')
    # Source Table 2 exact entries and boundary-invariance checks.
    assert step((Q(0), Q(0), Q(0))) == (Q(0), Q(0), Q(0))
    assert step((Q(1), Q(0), Q(0))) == (Q(1), Q(0), Q(0))
    for r in ((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)),
              (Q(0), Q(0), Q(1)), (Q(3, 5), Q(4, 5), Q(0))):
        assert norm_sq(step(r)) == 1
    for k in range(1, 7):
        r = (Q(0), Q(0), 1 - Q(1, 10 ** k))
        assert norm_sq(step(r)) > norm_sq(r)  # Purity need not decrease.
        assert 1 - norm_sq(step(r)) > 0
    print('BOUNDARIES: exact C0/C1, pure-sphere examples, near-pure transient purification')


def main():
    print('ARXIV: 1808.08151v1; all evidence uses exact rational/symbolic algebra')
    print('SYMPY_VERSION:', s.__version__)
    symbolic_checks()
    rational_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
