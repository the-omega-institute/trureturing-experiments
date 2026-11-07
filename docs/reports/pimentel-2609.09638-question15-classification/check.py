#!/usr/bin/env python3
"""Exact helper checks for Pimentel, arXiv:2609.09638v1, Question 15.

These checks do not prove the universal classification. Decimal displays are
readings of expressions already compared exactly, never inputs to assertions.
Requires Python 3 and SymPy.
"""
from math import gcd
import platform
import sympy as s


def exact_zero(expr):
    expr = s.trigsimp(s.simplify(expr))
    assert not expr.atoms(s.Float), expr
    assert expr == 0, expr


def source_six_arc():
    """Reproduce §5.1's six-arc curve, unit speed, closure, and length."""
    t = s.symbols("t", real=True)
    root3 = s.sqrt(3)
    circles = [s.Matrix([s.cos(t) / 2, s.sin(t) / 2, root3 / 2]),
               s.Matrix([s.cos(t), s.sin(t), 0]),
               s.Matrix([s.cos(t) / 2, s.sin(t) / 2, -root3 / 2])]
    Q = s.Matrix([[s.Rational(1, 2), 0, root3 / 2],
                  [0, 1, 0], [-root3 / 2, 0, s.Rational(1, 2)]])
    assert s.simplify(Q.T * Q) == s.eye(3)
    curves = [circles[0].subs(t, 2 * t),
              Q * circles[1].subs(t, t - 3 * s.pi / 2),
              circles[2].subs(t, 2 * t - 3 * s.pi),
              Q * circles[2].subs(t, 4 * s.pi - 2 * t),
              circles[1].subs(t, 7 * s.pi / 2 - t),
              Q * circles[0].subs(t, 7 * s.pi - 2 * t)]
    bounds = [0, s.pi / 2, 3 * s.pi / 2, 2 * s.pi,
              5 * s.pi / 2, 7 * s.pi / 2, 4 * s.pi]
    for i, c in enumerate(curves):
        exact_zero(c.dot(c) - 1)
        exact_zero(c.diff(t).dot(c.diff(t)) - 1)
        j = (i + 1) % 6
        end = c.subs(t, bounds[i + 1])
        start = curves[j].subs(t, bounds[j])
        assert all(s.simplify(a - b) == 0 for a, b in zip(end, start))
        tangent = c.diff(t).subs(t, bounds[i + 1])
        tangent_next = curves[j].diff(t).subs(t, bounds[j])
        assert all(s.simplify(a - b) == 0 for a, b in zip(tangent, tangent_next))
    exact_zero(sum(bounds[i + 1] - bounds[i] for i in range(6)) - 4 * s.pi)
    print("FIDELITY: source six arcs lie on S^2, join C^1, close, have unit speed and length 4*pi")


def endpoint_families():
    """Endpoints are labels j modulo 2n at signed colatitude (2j+1)π/(2n)."""
    cases = 0
    for n in range(2, 33):
        for k in range(1, n):
            remaining = set(range(2 * n))
            components = 0
            while remaining:
                components += 1
                pending = [next(iter(remaining))]
                while pending:
                    j = pending.pop()
                    if j not in remaining:
                        continue
                    remaining.remove(j)
                    # Upper fixed arcs and lower arcs rotated through 2kθ_n.
                    pending += [(-j - 1) % (2 * n), (2 * k - j - 1) % (2 * n)]
            assert components == gcd(n, k), (n, k, components)
            cases += 1
    print(f"FAMILY: {cases} endpoint graphs for 2<=n<=32 have exactly gcd(n,k) loops")


def lengths_and_means():
    """Length identity is checked in exact cyclotomic fields, not by decimals."""
    z = s.symbols("z")
    for n in range(2, 13):
        order = 4 * n
        numer = sum(z ** (2 * i + 1) - z ** ((-2 * i - 1) % order)
                    for i in range(n))
        # sin(theta)*sum_i sin((2i+1)*theta)=1, at z=exp(i*theta).
        residual = s.expand((z - z ** (order - 1)) * numer + 4)
        assert s.rem(residual, s.cyclotomic_poly(order, z), domain=s.QQ) == 0
    print("FAMILY: sum_i sin((2i+1)*theta_n)=1/sin(theta_n), exactly for 2<=n<=12")
    r = s.symbols("r", real=True)
    theta = s.symbols("theta", positive=True)
    primitive = r + s.cos(r) / s.sin(theta)
    exact_zero(s.diff(primitive, r) - (1 - s.sin(r) / s.sin(theta)))
    integral = primitive.subs(r, theta) - primitive.subs(r, 0)
    for n, expected_L, expected_J in [
        (2, 2 * s.pi * s.sqrt(2), s.pi / 4 - s.sqrt(2) + 1),
        (3, 4 * s.pi, s.pi / 6 + s.sqrt(3) - 2),
        (4, 2 * s.pi * s.sqrt(4 + 2 * s.sqrt(2)),
         s.pi / 8 + 1 + s.sqrt(2) - s.sqrt(4 + 2 * s.sqrt(2)))
    ]:
        th = s.pi / (2 * n)
        exact_zero(2 * s.pi / s.sin(th) - expected_L)
        exact_zero(integral.subs(theta, th) - expected_J)
        exact_zero(expected_J - (th - s.tan(th / 2)))
        # β's source tube law supplies this integral; no sampled quadrature.
        for k in range(1, n):
            if gcd(k, n) == 1:
                print(f"EXACT beta^{{{n},{k}}}: L={expected_L}; J={expected_J}")
        if n in (2, 3):
            print(f"DECIMAL DISPLAY ONLY n={n}: L={s.N(expected_L, 25)}, J(beta)={s.N(expected_J, 25)}")
    print("FIDELITY: source table for n=2,3,4 agrees exactly with beta tube-layer integral")


def rigidity_perturbation():
    """Check the non-binding lemma's local linear-algebra construction exactly."""
    p1 = s.Matrix([s.Rational(4, 5), 0, s.Rational(3, 5)])
    p2 = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    p3 = s.Matrix([-s.Rational(4, 5), 0, s.Rational(3, 5)])
    p4 = s.Matrix([0, -1, 0])
    u = s.Matrix([0, 0, 1])
    d = s.Rational(3, 5)
    M = s.Matrix.vstack(p1.T, p2.T, p3.T)
    assert M.det() != 0
    w = M.inv() * s.Matrix([1, -1, 1])
    eps = s.Rational(1, 100)
    v = u + eps * w
    margins = [s.simplify(v.dot(p) - d) for p in [p1, p2, p3, p4]]
    for p in [p1, p2, p3, p4]:
        exact_zero(p.dot(p) - 1)
    assert [s.sign(m) for m in margins] == [1, -1, 1, -1]
    radius_squared = 1 - d**2 / v.dot(v)
    assert 0 < radius_squared < s.Rational(3, 4)  # theta=pi/3
    print(f"RIGIDITY HELPER: det={M.det()}, w={list(w)}, strict alternating margins={margins}")
    print(f"RIGIDITY HELPER: perturbed circle radius^2={radius_squared} < sin(pi/3)^2=3/4")
    print("SCOPE: helper identities and finite family construction checked; no universal classification certified")


def main():
    print(f"Python {platform.python_version()}; SymPy {s.__version__}")
    source_six_arc()
    endpoint_families()
    lengths_and_means()
    rigidity_perturbation()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
