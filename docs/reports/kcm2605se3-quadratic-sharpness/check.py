#!/usr/bin/env python3
"""Exact source-fidelity and proposed se(3) sharpness checks; no numeric floats.

Source: Sooraj K.C and Vivek Mishra, arXiv:2605.31076v1.
These symbolic and finite checks support a preregistration, not a Lean proof.
"""
from fractions import Fraction
import sympy as sp


def equal(a, b):
    assert sp.simplify(sp.trigsimp(a - b)) == 0, (a, b)


t, r, s = sp.symbols("t r s", real=True)
R = sp.symbols("R", positive=True)

# Fidelity: Section 2.1 distinguishes the six-vector and matrix Frobenius norms.
w0, w1, w2, v0, v1, v2 = sp.symbols("w0 w1 w2 v0 v1 v2", real=True)
X = sp.Matrix([[0, -w2, w1, v0], [w2, 0, -w0, v1],
               [-w1, w0, 0, v2], [0, 0, 0, 0]])
frob_sq = sp.trace(X.T * X)
equal(frob_sq, 2 * (w0**2 + w1**2 + w2**2) + v0**2 + v1**2 + v2**2)
print("FIDELITY eq.(2): ||X||_mat^2 = 2|omega|^2 + |v|^2")

# Fidelity: Section 3.1 explicitly states these removable-singularity limits.
equal(sp.limit((1 - sp.cos(t)) / t**2, t, 0), sp.Rational(1, 2))
equal(sp.limit((t - sp.sin(t)) / t**3, t, 0), sp.Rational(1, 6))
print("FIDELITY Section 3.1: Rodrigues coefficient limits = 1/2, 1/6")

# Compute a fresh exact exponential on the actual rotation/translation family.
Xt = X.subs({w0: 0, w1: 0, w2: t, v0: r, v1: 0, v2: 0})
Et = Xt.exp()
q = sp.simplify(Et[1, 3])
equal(q, r * (1 - sp.cos(t)) / t)
equal(Et[0, 3], r * sp.sin(t) / t)
equal(Et[0, 0], sp.cos(t))
equal(Et[1, 0], sp.sin(t))
equal(sp.limit(q, t, 0), 0)
equal(sp.limit(sp.diff(q, t), t, 0), r / 2)
f = sp.cos(q)
fprime0 = sp.limit(sp.diff(f, t), t, 0)
q_jet = sp.series(q, t, 0, 4).removeO()
equal(q_jet, r*t/2 - r*t**3/24)
f_jet = sp.series(sp.cos(q_jet), t, 0, 3).removeO()
fsecond0 = sp.diff(f_jet, t, 2).subs(t, 0)
equal(fprime0, 0)
equal(fsecond0, -r**2 / 4)
equal(fsecond0.subs(r, R / 2), -R**2 / 16)
print("FAMILY: exp(X(t))_24 = r(1-cos(t))/t, smooth value 0 at t=0")
print("LOWER: f'(0)=0; f''(0)=-r^2/4=-R^2/16 at r=R/2")

# Global ambient derivatives: the only nonzero coordinate is G_24.
g = sp.symbols("g0:16", real=True)
F = sp.cos(g[7])  # row 2, column 4, in one-based row-major indexing.
grad = sp.Matrix([sp.diff(F, x) for x in g])
hess = sp.hessian(F, g)
assert all(grad[i] == 0 for i in range(16) if i != 7)
equal(grad[7], -sp.sin(g[7]))
assert all(hess[i, j] == 0 for i in range(16) for j in range(16)
           if (i, j) != (7, 7))
equal(hess[7, 7], -sp.cos(g[7]))
equal(1 - (grad.T * grad)[0], sp.cos(g[7])**2)
# H is rank one: H^T H = cos(g_24)^2 times the coordinate projection.
P = sp.zeros(16); P[7, 7] = 1
assert (hess.T * hess - sp.cos(g[7])**2 * P).applyfunc(sp.simplify) == sp.zeros(16)
equal(1 - sp.cos(g[7])**2, sp.sin(g[7])**2)
print("ADMISSIBILITY: 1-||grad F||_F^2=cos(g_24)^2 >= 0")
print("ADMISSIBILITY: H^T H=cos(g_24)^2 P_24; 1-cos(g_24)^2=sin(g_24)^2 >= 0")
print("ADMISSIBILITY: M1=M2=1 globally; F is a cosine of a linear coordinate")

# Exact simplex volumes in the Duhamel derivative bounds.
u, a = sp.symbols("u a", real=True)
first = sp.integrate(s, (s, 0, 1))
second = sp.integrate(2 * sp.integrate(sp.integrate(1, (u, 0, a)), (a, 0, s)), (s, 0, 1))
equal(first, sp.Rational(1, 2))
equal(second, sp.Rational(1, 3))
print("UPPER ingredients: integral s ds=1/2; two ordered simplex volumes integrate to 1/3")
base = R**2 / 4 + R / 3 + sp.Rational(3, 2)
upper = sp.Rational(25, 12) * R**2
equal(upper - base, (R - 1) * (11 * R + 9) / 6)
print("UPPER: 25R^2/12-(R^2/4+R/3+3/2)=(R-1)(11R+9)/6")

for radius in [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3),
               Fraction(10), Fraction(100)]:
    lower = radius**2 / 16
    hessian_bound = radius**2 / 4 + radius / 3 + Fraction(3, 2)
    upper_bound = Fraction(25, 12) * radius**2
    assert 0 < lower <= hessian_bound <= upper_bound
    # t=R/2 gives an interior point of the same six-dimensional ball.
    r0, t0 = radius / 2, radius / 2
    assert r0*r0 + t0*t0 < radius*radius
    print(f"SAMPLE R={radius}: lower={lower}; Hessian upper={hessian_bound}; quadratic upper={upper_bound}")

print("SCOPE: exact symbolic identities and rational samples; global derivative estimates and supremum passage require analytic proof")
print("ALL_CHECKS_OK")
