#!/usr/bin/env python3
"""Exact checks for the singular four-Gaussian sign identity.

All matrix entries and comparisons are exact SymPy expressions.  The three
configurations have ranks one, two, and three; their direct spherical-area
values are compared with the correlation formula.  The first check also
reproduces the source's Theorem 7.1 value at a right angle.
"""

from itertools import combinations, product
from fractions import Fraction
import sympy as sp


def sigma(x):
    return sp.Integer(1) if sp.ask(sp.Q.nonnegative(x)) is not False and x >= 0 else sp.Integer(-1)


def principal_minor_sum(R, k):
    n = R.rows
    if k == 0:
        return sp.Integer(1)
    total = sp.Integer(0)
    for I in combinations(range(n), k):
        total += R.extract(I, I).det()
    return sp.simplify(total)


def selector(R):
    e = [principal_minor_sum(R, k) for k in range(4)]
    rank = 3 if e[3] > 0 else (2 if e[2] > 0 else 1)
    Q = sp.zeros(4)
    for k in range(rank + 1):
        Q += (-1) ** k * e[rank - k] * (R ** k)
    cols = [j for j in range(4) if Q[j, j] > 0]
    assert cols, (R, Q)
    j = min(cols)
    a = Q[:, j]
    assert all(sp.simplify(x) == 0 for x in R * a)
    return rank, Q, a


def formula(R):
    rank, Q, a = selector(R)
    eps = [sigma(x) for x in a]
    P = sp.prod(eps)
    pair = sum(P * eps[i] * eps[j] * sp.asin(R[i, j])
               for i in range(4) for j in range(i + 1, 4))
    return rank, sp.simplify(-P - sp.Rational(2, 1) / sp.pi * pair)


def source_second_brightness_value():
    # Theorem 7.1 in the source: Omega(theta) at theta = pi/2.
    theta = sp.pi / 2
    omega = sp.Rational(2, 3) / sp.pi * (
        sp.sin(theta) + sp.cos(theta) * (sp.pi / 2 - theta)
    )
    assert sp.simplify(omega - sp.Rational(2, 3) / sp.pi) == 0
    return omega


def check_four_sign_identity():
    # On the 14 nonconstant sign patterns, product(t) = -1 - sum_{i<j} t_i t_j.
    checked = 0
    for t in product((-1, 1), repeat=4):
        if len(set(t)) == 1:
            continue
        lhs = sp.prod(t)
        rhs = -1 - sum(t[i] * t[j] for i in range(4) for j in range(i + 1, 4))
        assert lhs == rhs
        checked += 1
    assert checked == 14
    return checked


def gram(vectors):
    return sp.Matrix([[sp.simplify((vectors[i].T * vectors[j])[0]) for j in range(4)]
                      for i in range(4)])


def check_rank_one():
    e = sp.Matrix([1, 0, 0])
    R = gram([e, e, -e, -e])
    rank, value = formula(R)
    # All four great circles coincide; the sign product is identically +1.
    direct = sp.Integer(1)  # signed area 4*pi divided by sphere area 4*pi
    assert rank == 1 and sp.simplify(value - direct) == 0
    return rank, value, direct


def check_rank_two():
    e1 = sp.Matrix([1, 0, 0])
    e2 = sp.Matrix([0, 1, 0])
    q = sp.sqrt(2) / 2
    vectors = [e1, e2, sp.Matrix([q, q, 0]), sp.Matrix([q, -q, 0])]
    R = gram(vectors)
    rank, value = formula(R)
    # Four lines on the great circle have successive angular gaps pi/4.
    # The product signs on the eight arcs are +,-,+,-,+,-,+,-, so the
    # signed spherical area is 0 and the expectation is exactly 0.
    arc_signs = [1, -1, 1, -1, 1, -1, 1, -1]
    direct = sp.simplify(sum(arc_signs) * sp.pi / 4 / sp.pi)
    assert rank == 2 and sp.simplify(value - direct) == 0
    return rank, value, direct


def check_rank_three():
    e1 = sp.Matrix([1, 0, 0])
    e2 = sp.Matrix([0, 1, 0])
    e3 = sp.Matrix([0, 0, 1])
    q = sp.sqrt(3) / 3
    vectors = [e1, e2, e3, sp.Matrix([q, q, q])]
    R = gram(vectors)
    rank, value = formula(R)
    alpha = sp.asin(sp.sqrt(3) / 3)
    # In one one-negative octant the cut x+y+z=0 bounds a quadrilateral
    # with vertices e1, e2, (0,1,-1)/sqrt(2), (1,0,-1)/sqrt(2).  Triangulate
    # it and use the exact solid-angle formula
    # 2*atan(|det(a,b,c)|/(1+a.b+b.c+c.a)).  The two triangle areas are
    # pi/4 and 2*atan((3-sqrt(2))/7), whose tangent addition identity is
    # exactly the tangent of 2*asin(1/sqrt(3)); hence A=2*alpha.
    A0 = sp.Matrix([1, 0, 0])
    A1 = sp.Matrix([0, 1, 0])
    A2 = sp.Matrix([0, 1, -1]) / sp.sqrt(2)
    A3 = sp.Matrix([1, 0, -1]) / sp.sqrt(2)

    def spherical_triangle_area(a, b, c):
        det = sp.Abs(sp.det(sp.Matrix.hstack(a, b, c)))
        den = 1 + a.dot(b) + b.dot(c) + c.dot(a)
        return 2 * sp.atan(det / den)

    triangle_one = spherical_triangle_area(A0, A1, A2)
    triangle_two = spherical_triangle_area(A0, A2, A3)
    x = (3 - sp.sqrt(2)) / 7
    tan_sum = (1 + 2 * x / (1 - x**2)) / (1 - 2 * x / (1 - x**2))
    assert sp.simplify(triangle_one - sp.pi / 4) == 0
    assert sp.simplify(triangle_two - 2 * sp.atan(x)) == 0
    assert sp.simplify(tan_sum - 2 * sp.sqrt(2)) == 0
    # In one one-negative octant, the cut x+y+z=0 bounds a spherical
    # quadrilateral of area A=2*alpha.
    # Both sides are in (0, pi), so the exact tangent identity selects the
    # principal branch.  This is the only inverse-trigonometric branch choice.
    quadrilateral_area = 2 * alpha
    # There are six one/two-negative octants.  The two all-equal octants
    # contribute signed area pi, hence total signed area 4*pi-24*alpha.
    direct = sp.simplify((4 * sp.pi - 12 * quadrilateral_area) / (4 * sp.pi))
    assert rank == 3 and sp.simplify(value - direct) == 0
    return rank, value, direct


def main():
    print("source_theorem_7_1_omega_pi_over_2 =", source_second_brightness_value())
    print("finite_sign_patterns_checked =", check_four_sign_identity())
    for name, check in (("rank1", check_rank_one), ("rank2", check_rank_two),
                        ("rank3", check_rank_three)):
        rank, value, direct = check()
        print(name, "rank =", rank, "formula =", value, "direct_area =", direct)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
