#!/usr/bin/env python3
"""Exact algebra checks for Bialy's magnetic rotational-rigidity question.

This checks source identities and the proposed local geometric mechanism.
It does not prove the global invariant-circle or twist-map assertions.
"""

import sympy as sp


def zero(expr):
    reduced = sp.simplify(sp.expand(expr))
    assert reduced == 0, reduced


def dot(a, b):
    return (a.T * b)[0]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def norm2(a):
    return dot(a, a)


J = sp.Matrix([[0, -1], [1, 0]])


def source_fidelity():
    R, h, hp, rho, x, r = sp.symbols("R h hp rho x r", positive=True)
    # Section 4, the polynomial in cos(t) in the Jacobian formula and its endpoint lemma.
    f = ((h**2 + rho * h + hp**2) * x**2
         - (2 * R * h + rho * (R**2 + h**2) / R) * x
         + (R**2 - hp**2 + rho * h))
    zero(f.subs(x, 1) - (R - h)**2 * (R - rho) / R)
    zero(f.subs(x, -1) - (R + h)**2 * (R + rho) / R)
    fc = f.subs({h: r, hp: 0, rho: r})
    zero(fc - 2 * r**2 * (x - R/r) * (x - (R**2 + r**2)/(2*r*R)))
    zero(fc.subs(x, 1) - (R-r)**3 / R)
    print("FIDELITY: both source endpoint identities and the circular factorization OK")


def standard_twist_system():
    # Solve the actual two equations displayed in Section 3 before t'_beta.
    # a = delta-t, b = s; sin(mu) = sin(a-b) = sa*cb-ca*sb.
    R, rho, sa, ca, sb, cb = sp.symbols("R rho sa ca sb cb", real=True)
    tp, one_minus_sp = sp.symbols("tp one_minus_sp", real=True)
    equations = [sa * tp - sb * one_minus_sp,
                 R * ca * tp + rho - R * cb * one_minus_sp]
    solution = sp.solve(equations, [tp, one_minus_sp], dict=True)[0]
    expected = rho * sb / (R * (sa * cb - ca * sb))
    zero(solution[tp] - expected)
    for eq in equations:
        zero(eq.subs(solution))
    print("TWIST SYSTEM: t'_beta = rho_beta*sin(s)/(R*sin(mu)) OK")
    print("SOURCE CORRECTION: printed t'_beta has wrong sign and lacks 1/R; G_12 lacks cancellation of R")


def symbolic_radial_jump():
    c, s, x, y, u, v, R = sp.symbols("c s x y u v R", real=True)
    A = sp.Matrix([[c, -s], [s, c]])
    q, C = sp.Matrix([x, y]), sp.Matrix([u, v])
    # For an actual rotation, A^-1 = A^T. Successive centers are C and A^T C.
    jump = -J * (C - A.T * C) / R
    bisector = dot(C, A*q-q)
    zero(R * det(q, jump) - bisector)
    # Equal Larmor radii give the bisector equation, modulo c^2+s^2=1.
    difference = norm2(C-A*q) - norm2(C-q)
    zero(difference + 2*bisector - (c*c+s*s-1)*norm2(q))
    print("SYMBOLIC: R*det(q, velocity_jump) = C dot (Aq-q) OK")
    print("SYMBOLIC: equal Larmor radii imply the radial-jump identity for every rotation OK")


def exact_regular_orbits():
    count = 0
    # All nontrivial rotation steps at these orders, including nonprimitive steps.
    # All quantities are integers, rationals, or exact quadratic radicals.
    for k in (2, 3, 4, 6):
        for m in range(1, k):
            theta = 2 * sp.pi * sp.Rational(m, k)
            c, s = sp.cos(theta), sp.sin(theta)
            A = sp.Matrix([[c, -s], [s, c]])
            zero(norm2(A * sp.Matrix([1, 0])) - 1)
            for r in (1, 2):
                for R in range(r+1, r+4):
                    q = sp.Matrix([r, 0])
                    qnext = A*q
                    chord = qnext-q
                    d2 = sp.simplify(norm2(chord))
                    assert d2 > 0
                    offset = sp.sqrt(sp.Rational(R*R, 1)/d2 - sp.Rational(1, 4))
                    C = sp.simplify((q+qnext)/2 + offset*J*chord)
                    zero(norm2(C-q) - R*R)
                    zero(norm2(C-qnext) - R*R)
                    vout = sp.simplify(-J*(C-q)/R)
                    vin = sp.simplify(-J*(A.T*C-q)/R)
                    zero(norm2(vout)-1)
                    zero(norm2(vin)-1)
                    assert sp.simplify(dot(vout, q)) < 0  # leaves the circle inward
                    jump = sp.simplify(vout-vin)
                    zero(det(q, jump))
                    assert sp.simplify(norm2(jump)) > 0  # excludes grazing
                    n = -q/r
                    reflection = sp.simplify(vin-2*dot(n, vin)*n)
                    for z in reflection-vout:
                        zero(z)
                    # Reflection at the next vertex is its rotated copy.
                    incoming_next = sp.simplify(-J*(C-qnext)/R)
                    nnext = -qnext/r
                    for z in incoming_next-2*dot(nnext, incoming_next)*nnext-A*vout:
                        zero(z)
                    count += 1
    print(f"EXACT CASES: {count} circular orbit configurations; orders 2,3,4,6; every nonzero step; two table radii and three weak-field radii each OK")


def noncircular_sensitivity():
    # A noncircular ellipse has a nonradial normal at a generic exact point.
    # Thus the proposed radial-normal conclusion does not hold for arbitrary tables.
    q = sp.Matrix([sp.Rational(6, 5), sp.Rational(4, 5)])
    normal = sp.Matrix([q[0]/4, q[1]])  # ellipse x^2/4 + y^2 = 1
    zero(q[0]**2/4 + q[1]**2 - 1)
    assert det(q, normal) == sp.Rational(18, 25)
    print("SENSITIVITY: noncircular ellipse normal is nonradial (det = 18/25) OK")


if __name__ == "__main__":
    print(f"Python exact arithmetic; SymPy {sp.__version__}; no floating point")
    source_fidelity()
    standard_twist_system()
    symbolic_radial_jump()
    exact_regular_orbits()
    noncircular_sensitivity()
    print("ALL_CHECKS_OK")
