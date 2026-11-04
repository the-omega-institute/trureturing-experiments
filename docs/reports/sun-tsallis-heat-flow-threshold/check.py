#!/usr/bin/env python3
"""Exact algebra for Sun's heat-flow lower-endpoint question.

Requires Python 3 and SymPy. No floating point and no numerical quadrature.
These checks are evidence for a preregistration, not an analytic or Lean proof.
"""
from fractions import Fraction
import sympy as sp


def zero(expr):
    assert sp.simplify(expr) == 0, expr


# Fidelity: Sun arXiv:2601.12944v3, Remark 20, for the kernel with variance 2t.
d, q, t = sp.symbols('d q t', positive=True)
alpha = d * (1 - q) / 2
U = q ** (-d / 2) * (4 * sp.pi * t) ** alpha
S = (1 - U) / (q - 1)
source_second = d / (2 * t**2) * (d * (1 - q) / 2 - 1) * U
zero(sp.diff(S, t, 2) - source_second)
assert sp.simplify(source_second.subs({d: 3, q: sp.Rational(1, 4), t: 1})).is_positive
print('FIDELITY: Remark 20 Gaussian entropy second derivative agrees exactly')

# Algebraic square completion, conditional on the stated IBP identities.
D, A, B, delta = sp.symbols('D A B delta', real=True)
c = (d - 1)**2 / (d * (d + 2))
HH = (d - 1) * D / d
VV = (d - 1) * B / d
HV = B / 2 - (d + 2) * A / (2 * d)
Q = HH - 2 * d / (d - 1) * HV + d**2 / (d - 1)**2 * VV
zero(A + c * D - (d - 1) / (d + 2) * Q)
zero(D + delta * A - ((1 - delta * c) * D + delta * (d - 1) / (d + 2) * Q))
qd = 2 * (d - 1)**2 / (2 * d**2 + 1)
zero((1 - (2 / q - 1) * c) - (2 * d**2 + 1) * (q - qd) / (q * d * (d + 2)))
zero(qd - (1 - 2 / d) - (d + 2) / (d * (2 * d**2 + 1)))
print('IDENTITY: square coefficients and exact threshold factorization agree')

s = sp.symbols('s', nonnegative=True)
for dim in range(2, 13):
    m = dim - 1
    u = (1 + s)**(-m)
    us = sp.diff(u, s)
    uss = sp.diff(u, s, 2)
    lap = 4 * s * uss + 2 * dim * us
    zero(lap - 2 * dim * m * (s - 1) * (1 + s)**(-dim - 1))
    # For radial u(s), trace-free Hessian and V both are multiples of
    # x tensor x - (s/d) Id. Verify the coefficient, for every radius s.
    zero(4 * uss - sp.Rational(dim, m) * 4 * us**2 / u)
    # The angular factor |S^(d-1)|/2 is common and strictly positive.
    N = 2 * dim + 2
    def radial_integral(poly):
        ans = 0
        for (power,), coefficient in sp.Poly(poly, s).terms():
            a = sp.Rational(dim, 2) + power
            assert 0 < a < N
            ans += coefficient * sp.gamma(a) * sp.gamma(N - a) / sp.gamma(N)
        return sp.simplify(ans)
    exact_D = 4 * dim**2 * m**2 * radial_integral((s - 1)**2)
    exact_A = 8 * dim * m**3 * radial_integral(s * (s - 1))
    cd = Fraction(m*m, dim*(dim+2))
    endpoint = Fraction(2*m*m, 2*dim*dim+1)
    assert exact_D.is_positive
    zero(exact_A / exact_D + sp.Rational(cd.numerator, cd.denominator))
    for ratio in (Fraction(1, 2), Fraction(9, 10), Fraction(99, 100)):
        index = ratio * endpoint
        exponent = 2 * Fraction(m, 1) / index
        # Normalizability, finite initial q-entropy, and strict initial S''.
        assert exponent > Fraction(dim, 2)
        assert 2*m > Fraction(dim, 2)
        second_per_D = -8 / index * (1 - (2 / index - 1) * cd)
        assert second_per_D > 0
    zero(1 - (2 / sp.Rational(endpoint.numerator, endpoint.denominator) - 1) * sp.Rational(cd.numerator, cd.denominator))
    print(f'd={dim}: q_d={endpoint}; A/D={-cd}; D/angular={exact_D}; three strict below-endpoint signs OK')

print('LIMIT: heat-convolution estimates, C2 passage, and positive-time persistence require analytic proof')
print('ALL_CHECKS_OK')
