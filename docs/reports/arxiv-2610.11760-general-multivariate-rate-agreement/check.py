#!/usr/bin/env python3
"""Exact symbolic certificates for the Scricciolo Section 8 candidate.

Source: arXiv:2610.11760v1, Definition 2.1 and Section 8.
Uses SymPy exact algebra and symbolic integration, without floating point.
The statistical and measure-theoretic deductions are explained in README.md;
this program checks their explicit analytic and algebraic inputs, not Lean proofs.
"""

import sympy as s


def equal(a, b):
    assert s.simplify(a - b) == 0, (a, b)


R = s.Rational
x, t = s.symbols("x t", real=True)
a = s.symbols("a", nonnegative=True)
n = s.symbols("n", positive=True, integer=True)
theta = s.symbols("theta", real=True)

# Fidelity: the source explicitly gives the standard Laplace CF and reciprocal.
cf = s.integrate(s.exp(-x) * s.exp(s.I * t * x) / 2, (x, 0, s.oo), conds="none")
cf += s.integrate(s.exp(x) * s.exp(s.I * t * x) / 2, (x, -s.oo, 0), conds="none")
# The half-line integral convergence conditions are Re(1 +/- i*t)>0.
equal(s.re(1 + s.I*t), 1)
equal(s.re(1 - s.I*t), 1)
equal(cf, 1 / (1 + t**2))
equal(s.diff(1 + t**2, t), 2 * t)
equal(s.diff(1 + t**2, t, 2), 2)
equal(s.integrate(x**2 * s.exp(-x), (x, 0, s.oo), conds="none"), 2)
# Clearing positive denominators certifies c0=1, c1=2, beta=2.
equal((1 + a)**2 - (1 + a**2), 2 * a)
equal(2 * (1 + a**2) - (1 + a)**2, (a - 1)**2)
equal(2 * (1 + a)**2 - 2 * a, 2 * (1 + a + a**2))
equal(2 * (1 + a)**2 - 2, 2 * a * (2 + a))
E = 1 / (1 + 2 * (R(2) - R(1, 2)) / (R(1) + 1))
equal(E, R(2, 5))
equal(E / 2, R(1, 5))
print("FIDELITY: Laplace CF=1/(1+t^2), reciprocal derivatives=(2t,2), variance=2")
print("FIDELITY: source exponent at d=2, alpha=(1,1), beta=(2,2) is 1/5")

# Correlated Gaussian family: X_theta=(theta+U, U/2+sqrt(3)V/2).
Sigma = s.Matrix([[1, R(1, 2)], [R(1, 2), 1]])
B = s.Matrix([[1, 0], [R(1, 2), s.sqrt(3) / 2]])
assert B * B.T == Sigma
assert Sigma.eigenvals() == {R(1, 2): 1, R(3, 2): 1}
assert Sigma.det() == R(3, 4)
assert Sigma.inv() == s.Matrix([[R(4, 3), -R(2, 3)], [-R(2, 3), R(4, 3)]])
assert R(4, 3) < 4  # L1 derivative <= sqrt(4/3) < 2.
u, v = s.symbols("u v", real=True)


def normal_moment(k):
    return s.Integer(0) if k % 2 else s.Integer(1) if k == 0 else s.factorial2(k - 1)


def gaussian_expectation(poly):
    out = s.Integer(0)
    for (i, j), coefficient in s.Poly(s.expand(poly), u, v).terms():
        out += coefficient * normal_moment(i) * normal_moment(j)
    return s.expand(out)


X1, X2 = theta + u, (u + s.sqrt(3) * v) / 2
equal(gaussian_expectation(X1), theta)
equal(gaussian_expectation(X2), 0)
equal(gaussian_expectation(X1 * X2), R(1, 2))
equal(gaussian_expectation((1 + X1**2) * (1 + X2**2)), R(9, 2) + 2 * theta**2)
equal(R(13, 2) - (R(9, 2) + 2 * theta**2), 2 * (1 - theta) * (1 + theta))
print("FAMILY: Sigma eigenvalues={1/2,3/2}, covariance=1/2; dependent and nondegenerate")
print("CLASS: derivative squared bounds=4/3<4; N_1^2((1,1),2)")
print("MOMENT: E[(1+X1^2)(1+X2^2)]=9/2+2 theta^2<=13/2 for 0<=theta<=1")

# A positive rectangle gives a strip certificate without numerical CDF evaluation.
# For x in [0,1/2], y>=1/4, the event Z1 in [-1/4,0], Z2 in [0,1/4]
# lies in {x-1 < Z1 <= x, Z2 <= y}, so F0-F1 is bounded below by c.
z1, z2 = s.symbols("z1 z2", real=True)
q = z1**2 - z1*z2 + z2**2
assert max(q.subs({z1: i, z2: j}) for i in [-R(1, 4), 0]
           for j in [0, R(1, 4)]) == R(3, 16)
equal(q.subs({z1: -u/4, z2: v/4}), (u**2 + u*v + v**2) / 16)
equal(3 - (u**2 + u*v + v**2), (1-u**2) + (1-u*v) + (1-v**2))
assert -R(1, 4) > R(1, 2) - 1
c = s.exp(-R(1, 8)) / (16 * s.pi * s.sqrt(3))
assert c.is_positive
T = s.symbols("T", positive=True)
strip_bound = c * R(1, 2) * T
equal(s.diff(strip_bound, T), c / 2)
assert s.limit(strip_bound, T, s.oo) == s.oo
print("CDF: F0(x,y)-F1(x,y)>=exp(-1/8)/(16*pi*sqrt(3)) on [0,1/2] x [1/4,infinity)")
print("CDF: integral over strips of height T >= exp(-1/8)*T/(32*pi*sqrt(3)) -> infinity")
print("CDF RISK: equivalent positive finite-sample densities + infinite separation => R_CDF(n)=infinity")

# The sample-mean upper bound and a two-point lower bound have the same order.
delta = 1 / (2 * s.sqrt(n))
shift = s.Matrix([delta, 0])
KL = n * (shift.T * Sigma.inv() * shift)[0] / 2
equal(KL, R(1, 6))
equal(KL / 2, R(1, 12))
assert R(1, 12) < R(1, 4)  # Pinsker TV <= sqrt(1/12) < 1/2.
equal(delta / 8, 1 / (16 * s.sqrt(n)))
equal((1 + 2) / n, 3 / n)
lower = 1 / (16 * s.sqrt(n))
upper = s.sqrt(3 / n)
equal(upper / lower, 16 * s.sqrt(3))
print("W1: translation coupling + first mean => W1(mu_theta,mu_eta)=abs(theta-eta)")
print("W1 LOWER: delta=1/(2*sqrt(n)); observation KL<=1/6; TV<1/2; R_W1(n)>=1/(16*sqrt(n))")
print("W1 UPPER: clipped first-coordinate sample mean gives R_W1(n)<=sqrt(3/n)")
print("SCOPE: exact symbolic certificates and a mathematical outline; no Lean settlement proof")
print("ALL_CHECKS_OK")
