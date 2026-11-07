#!/usr/bin/env python3
"""Exact checks for Wu's two FORMULA conjectures in OEIS A005908.

Finite checks support a proposed all-index argument; they do not prove it.
Run with Python 3 and SymPy. No network access or floating-point arithmetic.
"""

from hashlib import sha256
from math import isqrt

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero(expression, message):
    require(sp.simplify(expression) == 0, message)


def floor_phi(m):
    """floor((1+sqrt(5))*m/2), evaluated exactly for m >= 0."""
    require(m >= 0, "floor_phi domain")
    return (m + isqrt(5 * m * m)) // 2


def coefficients(numerator, denominator, count):
    """Integer coefficients of N/D when D(0)=1, by coefficient recursion."""
    require(denominator[0] == 1, "unit constant coefficient")
    result = []
    for n in range(count):
        result.append(
            (numerator[n] if n < len(numerator) else 0)
            - sum(denominator[j] * result[n - j]
                  for j in range(1, min(n + 1, len(denominator))))
        )
    return result


def main():
    # %S/%T/%U of https://oeis.org/A005908/internal, revision 27.
    source_terms = [
        0, 1, 1, 2, 4, 9, 20, 46, 106, 245, 567, 1313, 3041, 7044,
        16317, 37798, 87559, 202831, 469860, 1088436, 2521375, 5840796,
        13530276, 31343052, 72606569, 168194019, 389623535, 902567761,
        2090809436, 4843386045,
    ]
    count = 4000
    a = [0, 1]
    while len(a) < count:
        a.append(floor_phi(a[-1]) + floor_phi(a[-2]))
    require(a[:len(source_terms)] == source_terms, "OEIS source fidelity")
    print("FIDELITY_OK: all 30 displayed OEIS terms, indices 0..29")

    # Independent floor criterion after removing the radical:
    # 2*k-m <= sqrt(5)*m < 2*k+2-m. Signs are checked before squaring.
    for m in range(10001):
        k = floor_phi(m)
        left, right = 2 * k - m, 2 * k + 2 - m
        require(left <= 0 or left * left <= 5 * m * m, "floor lower bound")
        require(right > 0 and 5 * m * m < right * right, "floor upper bound")
    print("FLOOR_OK: exact squared inequalities for m=0..10000")

    d_coeffs = [1, -2, -1, 0, 1, 1]
    u = coefficients([0, 1, -1, -1, -1], d_coeffs, count)
    v = coefficients([0, 1, -1, 0, -1], d_coeffs, count)
    require(u == a, "floor recursion versus rational coefficients")
    for n in range(count):
        require(v[n] == floor_phi(u[n]), f"companion floor at {n}")
        require(sum(d_coeffs[j] * a[n - j]
                    for j in range(min(n + 1, len(d_coeffs))))
                == ([0, 1, -1, -1, -1][n] if n < 5 else 0),
                f"generating-function coefficient at {n}")
    for n in range(count - 2):
        require(u[n + 2] == v[n + 1] + v[n], f"integer companion at {n}")
    for n in range(5, count):
        require(a[n] == 2 * a[n - 1] + a[n - 2] - a[n - 4] - a[n - 5],
                f"conjectured recurrence at {n}")
    print("FINITE_IDENTITY_OK: recurrence at every n=5..3999")
    print("FINITE_GF_OK: D*A=N at every degree 0..3999")
    print("FINITE_COMPANION_OK: v_n=floor(phi*u_n), n=0..3999")

    x, q, q_next = sp.symbols("x q q_next", real=True)
    phi = (1 + sp.sqrt(5)) / 2
    s = phi - 1
    c = 1 / sp.sqrt(5)
    d = 1 - 2*x - x**2 + x**4 + x**5
    u_gf = x * (1 - x - x**2 - x**3) / d
    v_gf = x * (1 - x - x**3) / d
    zero(s**2 + s - 1, "quadratic identity")
    zero(c - 1 / (1 + 2*s), "centering constant")
    zero(c**2 - sp.Rational(1, 5), "squared centering constant")
    zero(d - (1-x)*(1-phi*x-phi*x**2)*(1+s*x+s*x**2), "factorization")
    zero((1+x)*v_gf - (u_gf/x - 1), "companion rational identity")
    r_gf = phi*u_gf - v_gf
    zero(r_gf - s*x*(1+s*x)/((1-x)*(1+s*x+s*x**2)), "error cancellation")
    q_gf = r_gf - c/(1-x)
    zero(q_gf - (-c + (s-phi*c)*x)/(1+s*x+s*x**2), "centered error GF")
    print("SYMBOLIC_GF_OK: factorization, companion, error and centered-error identities")

    def energy(z, w):
        return w*w + s*z*w + s*z*z

    k = s - s*s/4
    e0 = energy(-c, s-c)
    zero(energy(q_next, -s*(q_next+q)) - s*energy(q, q_next), "energy contraction")
    zero(energy(q, q_next) - ((q_next+s*q/2)**2 + k*q*q), "square completion")
    zero(e0 - s*(s-c), "initial energy form 1")
    zero(e0 - (3-4*s)/5, "initial energy form 2")
    zero(k*c*c - s*e0 - s*(15*s-8)/20, "strict margin identity")
    zero((q_gf*(1+s*x+s*x**2)).subs(x, 0) + c, "q0")
    zero((s-c) + s*(-c) - (s-phi*c), "q1")
    # All comparisons below are exact comparisons of rational numbers or real radicals.
    require(sp.Rational(31, 15)**2 < 5, "sqrt(5)>31/15")
    require(5 < sp.Rational(5, 2)**2, "sqrt(5)<5/2")
    require(4 < 5 < 9, "2<sqrt(5)<3")
    require(s > sp.Rational(8, 15), "s>8/15")
    require(0 < s < sp.Rational(3, 4) < 1, "s domain")
    require(0 < c and 2*c < 1, "floor strip upper endpoint")
    require(e0 > 0 and k > 0 and k*c*c-s*e0 > 0, "positive strict energy margin")
    print("SYMBOLIC_ENERGY_OK: E_next=s*E, initial energy, square completion")
    print("EXACT_BOUNDS_OK: 8/15<s<3/4, 0<c, 2*c<1, K*c^2-s*E0=s*(15*s-8)/20>0")
    print("DATA_SHA256:", sha256((",".join(map(str, a)) + "\n").encode()).hexdigest())
    print("SCOPE: finite checks and exact algebra; no all-index Lean proof or settlement claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
