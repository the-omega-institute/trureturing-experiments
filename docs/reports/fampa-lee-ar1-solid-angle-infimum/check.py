#!/usr/bin/env python3
"""Exact checks for Fampa--Lee, arXiv:2610.12665v1, §8.1.

All mathematical inputs are SymPy integers/rationals.  The decimal readings are
only high-precision displays of exact SymPy expressions; no binary floating
point is used to decide a check.
"""

from fractions import Fraction
from itertools import combinations

import sympy as sp


def cross(a, b):
    return sp.Matrix(a).cross(sp.Matrix(b))


def dot(a, b):
    return sp.Matrix(a).dot(sp.Matrix(b))


def trihedral_angle(u, v, w):
    """Van Oosterom--Strackee formula, with exact radical arguments."""
    u, v, w = map(sp.Matrix, (u, v, w))
    nu, nv, nw = (sp.sqrt(x.dot(x)) for x in (u, v, w))
    numerator = sp.Abs(sp.det(sp.Matrix.hstack(u, v, w)))
    denominator = nu * nv * nw + u.dot(v) * nw + v.dot(w) * nu + w.dot(u) * nv
    return 2 * sp.atan(sp.simplify(numerator / denominator))


def exact_hull(points):
    """Monotone-chain hull of 2D SymPy-rational points."""
    pts = sorted(points, key=lambda z: (z[0], z[1]))

    def orient(a, b, c):
        return sp.simplify((b[0] - a[0]) * (c[1] - a[1]) -
                           (b[1] - a[1]) * (c[0] - a[0]))

    lower = []
    for p in pts:
        while len(lower) >= 2 and orient(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and orient(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def cone_angle(rays):
    """Exact spherical area by an exact cross-section hull and triangulation."""
    # q=(1,1,1) is positive on every ray used below for the checked samples.
    normalized = []
    for r in rays:
        s = sp.simplify(sum(r))
        if not sp.ask(sp.Q.positive(s)):
            raise AssertionError(f"cross-section functional is not positive: {r}")
        p = (sp.cancel(r[0] / s), sp.cancel(r[1] / s))
        normalized.append((p[0], p[1], tuple(r)))
    hull = exact_hull([(p[0], p[1]) for p in normalized])
    by_point = {(p[0], p[1]): p[2] for p in normalized}
    ordered = [by_point[p] for p in hull]
    if len(ordered) < 3:
        raise AssertionError("cone cross-section is not full-dimensional")
    area = sp.Integer(0)
    for i in range(1, len(ordered) - 1):
        area += trihedral_angle(ordered[0], ordered[i], ordered[i + 1])
    return sp.simplify(area)


def ar1_matrix(rho):
    return sp.Matrix([[1, -rho, 0], [-rho, 1 + rho**2, -rho],
                      [0, -rho, 1]])


def sigma_denominator(rho):
    r13 = 3 * rho**2 / (rho**4 + rho**2 + 1)
    r12 = rho * (rho**2 + 2) / sp.sqrt(2 * rho**6 + 3 * rho**4 + 3 * rho**2 + 1)
    return sp.Rational(1, 8) + (sp.asin(r13) + 2 * sp.asin(r12)) / (4 * sp.pi)


def predicted_limit(rho):
    return (sp.Rational(1, 8) + sp.atan(rho**2) / (4 * sp.pi)) / sigma_denominator(rho)


def generators(rho, eps):
    A = ar1_matrix(rho)
    xi = sp.Matrix([1 + rho**2 + eps, rho + eps, eps])
    b = sp.simplify(A * xi)
    L = []
    for i in range(3):
        vals = [xi[k] / rho ** abs(k - i) for k in range(3)]
        # Inputs used by the samples are exact rationals, so this comparison
        # is an exact SymPy comparison (no decimal ordering is involved).
        L.append(sp.simplify(-(1 - rho**2) * min(vals)))
    e = [sp.eye(3).col(i) for i in range(3)]
    v = [sp.simplify(b + L[i] * e[i]) for i in range(3)]
    return xi, b, L, e + v


def exact_structure_checks():
    rho, eps = sp.symbols("rho eps", positive=True)
    A = ar1_matrix(rho)
    Sigma = sp.Matrix([[rho ** abs(i - j) for j in range(3)] for i in range(3)])
    H = Sigma / (1 - rho**2)
    assert sp.simplify(A * H - sp.eye(3)) == sp.zeros(3)

    b_star = sp.Matrix([1, 0, -rho**2])
    assert (Sigma * b_star)[2] == 0
    assert sp.simplify((Sigma * b_star)[0] - (1 - rho**4)) == 0
    assert sp.simplify((Sigma * b_star)[1] - rho * (1 - rho**2)) == 0

    # Proposition 8.16's trihedral part: e1,e2,b* has tan(Omega/2)
    # equal to rho^2/(sqrt(1+rho^4)+1), hence Omega=atan(rho^2).
    u, v, w = sp.eye(3).col(0), sp.eye(3).col(1), b_star
    tangent = sp.simplify(abs(sp.det(sp.Matrix.hstack(u, v, w))) /
                          (sp.sqrt(u.dot(u)) * sp.sqrt(v.dot(v)) * sp.sqrt(w.dot(w)) +
                           u.dot(v) * sp.sqrt(w.dot(w)) + v.dot(w) * sp.sqrt(u.dot(u)) +
                           w.dot(u) * sp.sqrt(v.dot(v))))
    expected_tangent = rho**2 / (sp.sqrt(1 + rho**4) + 1)
    assert sp.simplify(tangent - expected_tangent) == 0
    assert sp.simplify(expected_tangent - (sp.sqrt(1 + rho**4) - 1) / rho**2) == 0

    # The denominator in Proposition 8.15 is the Gaussian orthant formula
    # applied to Sigma^2; verify its two correlation entries exactly.
    S2 = Sigma * Sigma
    assert sp.simplify(S2[0, 2] / sp.sqrt(S2[0, 0] * S2[2, 2]) -
                       3 * rho**2 / (rho**4 + rho**2 + 1)) == 0
    assert sp.simplify(S2[0, 1] / sp.sqrt(S2[0, 0] * S2[1, 1]) -
                       rho * (rho**2 + 2) /
                       sp.sqrt(2 * rho**6 + 3 * rho**4 + 3 * rho**2 + 1)) == 0

    # Exact boundary sequence from the proof: b_epsilon=A*xi_epsilon and
    # H*b_epsilon=xi_epsilon, with strictly positive xi_epsilon.
    xi = sp.Matrix([1 + rho**2 + eps, rho + eps, eps])
    b = sp.simplify(A * xi)
    assert sp.simplify(H * b - xi) == sp.zeros(3, 1)
    assert all(sp.simplify(xi[i] - (0 if i == 2 else 0)) != 0 for i in range(3))
    return True


def sample_readings():
    print("sample readings (exact rationals evaluated to 40 decimal digits):")
    for rho in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(4, 5)):
        limit = predicted_limit(rho)
        print(f"rho={rho}: predicted_limit={sp.N(limit, 40)}")
        for eps in (sp.Rational(1, 100), sp.Rational(1, 1000), sp.Rational(1, 10000)):
            xi, b, L, rays = generators(rho, eps)
            angle = cone_angle(rays)
            ratio = sp.simplify((angle / (4 * sp.pi)) / sigma_denominator(rho))
            if not sp.ask(sp.Q.positive(ratio - limit)):
                raise AssertionError("sample ratio did not exceed the proposed infimum")
            print(f"  eps={eps}: ratio={sp.N(ratio, 40)}")
            # Every native generator tends to b*.  For this rational sample
            # range the exact minimum in each L_i is the xi_3 term, so the
            # symbolic branch is written explicitly for the limit check.
            evar = sp.symbols("evar", positive=True)
            xi_sym = sp.Matrix([1 + rho**2 + evar, rho + evar, evar])
            b_sym = ar1_matrix(rho) * xi_sym
            L_sym = [-(1 - rho**2) * evar / rho**2,
                     -(1 - rho**2) * evar / rho,
                     -(1 - rho**2) * evar]
            b_star = sp.Matrix([1, 0, -rho**2])
            rays_sym = [b_sym + L_sym[i] * sp.eye(3).col(i) for i in range(3)]
            for r in rays_sym:
                assert sp.simplify((r - b_star).subs(evar, 0)) == sp.zeros(3, 1)


def main():
    exact_structure_checks()
    sample_readings()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
