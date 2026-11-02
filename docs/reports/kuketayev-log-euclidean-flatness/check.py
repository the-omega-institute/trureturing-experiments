#!/usr/bin/env python3
"""Exact algebra for the covariance functional in arXiv:2608.25288v1.

Requires Python 3 and SymPy. No floating point, simulation, or Lean proof.
The finite noise has the first four moments of N(0,t); it is not Gaussian.
"""

from fractions import Fraction
from itertools import combinations
from math import comb

import sympy as s


def equal(actual, expected):
    difference = s.factor(s.cancel(actual - expected))
    assert difference == 0, (actual, expected, difference)


U, V, r, q, h = s.symbols("U V r q h", real=True)
t = s.symbols("t", positive=True)
N = t**2 + 5*t + 4
scale = t + s.Rational(5, 2)


def reduce_q(expr):
    num, den = s.fraction(s.cancel(expr))
    assert not den.has(q)
    return s.cancel(s.rem(num, q**2 - (1-r**2), q) / den)


def expectation(poly, moments):
    return s.factor(sum(coefficient * moments[(a, b)]
                        for (a, b), coefficient in s.Poly(poly, U, V).terms()))


# Fidelity: Section 5.2 states V_C=1 for every Gaussian bivariate |r|<1.
# Compute it from Wick fourth moments in the observed X,Y coordinates.
X, Y = s.symbols("X Y", real=True)
gaussian_fourths = {(4, 0): 3, (3, 1): 3*r, (2, 2): 1+2*r**2,
                    (1, 3): 3*r, (0, 4): 3}
gaussian_if = (X*Y-r*(X**2+Y**2)/2)/(1-r**2)
gaussian_variance = sum(c*gaussian_fourths[(a, b)]
                        for (a, b), c in s.Poly(s.expand(gaussian_if**2), X, Y).terms())
equal(gaussian_variance, 1)
print("FIDELITY Section 5.2: Gaussian bivariate V(r) = 1 identically")


# Independent symbolic influence-function calculation, with q^2=1-r^2.
source_if = (U*(r*U+q*V)-r*(U**2+(r*U+q*V)**2)/2)/(1-r**2)
latent_if = r*(U**2-V**2)/2 + q*U*V
equal(reduce_q(source_if-latent_if), 0)
means = {(2, 0): 1, (0, 2): 1, (1, 1): 0}
equal(expectation(latent_if, means), 0)
m40, m04, m31, m13, m22 = s.symbols("m40 m04 m31 m13 m22")
generic = {(4, 0): m40, (0, 4): m04, (3, 1): m31,
           (1, 3): m13, (2, 2): m22}
generic_variance = reduce_q(expectation(s.expand(latent_if**2), generic))
formula = m22 + r**2*(m40+m04-6*m22)/4 + r*q*(m31-m13)
equal(generic_variance, formula)
print("IDENTITY v(r) = m22 + r^2*(m40+m04-6*m22)/4 + r*q*(m31-m13)")


# The finite latent law: B uniform on four integer points; independent
# G_i in {0,+sqrt(3t),-sqrt(3t)} with probabilities {2/3,1/6,1/6};
# R^2 in {1+d,1-d} equiprobably, d^2=9/(4N); (U,V)=R(B+G)/sqrt(scale).
B = ((2, 1), (-2, -1), (1, -2), (-1, 2))
noise_probabilities = (Fraction(2, 3), Fraction(1, 6), Fraction(1, 6))
assert sum(noise_probabilities) == 1
assert sum(Fraction(1, 8)*a*b
           for _ in B for _ in range(2)
           for a in noise_probabilities for b in noise_probabilities) == 1
equal(scale**2-N, s.Rational(9, 4))
equal(1+s.Rational(9, 4)/N, scale**2/N)
assert s.Poly(4*N-9, t).all_coeffs() == [4, 20, 7]  # d^2<1 for t>0


def noise_moment(k):
    if k == 0:
        return s.Integer(1)
    return s.Integer(0) if k % 2 else (3*t)**(k//2)/3


def unscaled_moment(a, b):
    return s.expand(sum(s.Rational(1, 4) * sum(
        comb(a, i)*comb(b, j) * b1**(a-i)*b2**(b-j)
        * noise_moment(i)*noise_moment(j)
        for i in range(a+1) for j in range(b+1)) for b1, b2 in B))


def latent_moment(a, b):
    degree = a+b
    if degree % 2:
        equal(unscaled_moment(a, b), 0)
        return s.Integer(0)
    # E R^(2k) = sum over even j of binom(k,j) d^j.
    k = degree//2
    radial = sum(comb(k, j)*(s.Rational(9, 4)/N)**(j//2)
                 for j in range(0, k+1, 2))
    return s.factor(unscaled_moment(a, b)*radial/scale**k)


moments = {(a, b): latent_moment(a, b)
           for a in range(5) for b in range(5-a)}
expected = {(0, 0): 1, (1, 0): 0, (0, 1): 0,
            (2, 0): 1, (0, 2): 1, (1, 1): 0,
            (2, 2): 1, (3, 1): 3/N, (1, 3): -3/N,
            (4, 0): 3-s.Rational(7, 2)/N, (0, 4): 3-s.Rational(7, 2)/N}
for indices, value in expected.items():
    equal(moments[indices], value)
print("ATOMIC 72 labels: total mass 1; centered; covariance I; all fourth moments exact")
print("MOMENTS m22=1; m31=3/N; m13=-3/N; m40=m04=3-7/(2N); N=t^2+5t+4")

atomic_variance = reduce_q(expectation(s.expand(latent_if**2), moments))
equal(atomic_variance, 1-s.Rational(7, 4)*r**2/N+6*r*q/N)
equal(atomic_variance.subs(r, 0).subs(q, 1), 1)
algebraic_variance = atomic_variance.subs(q, s.sqrt(1-r**2))
equal(s.diff(algebraic_variance, r).subs(r, 0), 6/N)
equal((6/N).subs(t, s.Rational(1, 2)), s.Rational(8, 9))
print("SETTLEMENT v_t(r)=1-7r^2/(4N)+6r*sqrt(1-r^2)/N; v_t'(0)=6/N>0")

# A rational parameter with rho'(0)=1 preserves the derivative 8/9.
rho = h/(1+h**2/4)
q_h = (1-h**2/4)/(1+h**2/4)
equal(rho**2+q_h**2, 1)
equal(s.diff(rho, h).subs(h, 0), 1)
rational_variance = s.factor(atomic_variance.subs({r: rho, q: q_h, t: s.Rational(1, 2)}))
equal(s.diff(rational_variance, h).subs(h, 0), s.Rational(8, 9))
print("RATIONAL normalized v(h) =", rational_variance)
print("RATIONAL v(0)=1; v'(0)=8/9; |h|<2; rho(h)=h/(1+h^2/4)")

for rv, qv, ev in [(Fraction(0), Fraction(1), Fraction(1)),
                   (Fraction(3, 5), Fraction(4, 5), Fraction(4, 3)),
                   (Fraction(-3, 5), Fraction(4, 5), Fraction(12, 25)),
                   (Fraction(5, 13), Fraction(12, 13), Fraction(5828, 4563))]:
    equal(rv*rv+qv*qv, 1)
    value = atomic_variance.subs({r: s.Rational(rv), q: s.Rational(qv), t: s.Rational(1, 2)})
    equal(value, s.Rational(ev))
    print(f"EXACT r={rv}; q={qv}; v={s.factor(value)}")


# Append independent symmetric ±1 coordinates. At I, the complete Jacobian
# sends XX^T to its off-diagonal products. Check their covariance, not only
# their diagonal variances, in dimensions 2 through 8.
def product_moment(p, indices):
    powers = [indices.count(i) for i in range(p)]
    return moments[powers[0], powers[1]] * (0 if any(k % 2 for k in powers[2:]) else 1)


for p in range(2, 9):
    pairs = list(combinations(range(p), 2))
    for i, pair in enumerate(pairs):
        equal(product_moment(p, pair), 0)
        for j, other in enumerate(pairs):
            equal(product_moment(p, pair+other), int(i == j))
    print(f"BASELINE p={p}: complete {len(pairs)}x{len(pairs)} off-diagonal covariance = I")

# Exact block isolation in the source's resolvent formula for D log.
a = s.symbols("a", nonnegative=True)
for p in (2, 3, 4):
    C = s.eye(p)
    C[0, 1] = C[1, 0] = r
    resolvent = (C+a*s.eye(p)).inv()
    H = s.Matrix(p, p, lambda i, j: s.Symbol(f"H{min(i,j)}{max(i,j)}"))
    full = (resolvent*H*resolvent)[0, 1]
    block = (resolvent[:2, :2]*H[:2, :2]*resolvent[:2, :2])[0, 1]
    equal(full, block)
print("BLOCK leading log derivative independent of appended coordinates (symbolic p=2,3,4)")
print("ALL_CHECKS_OK")
