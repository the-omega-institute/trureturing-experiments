#!/usr/bin/env python3
"""Exact checks for the Section 4.2 convergence conjecture in arXiv:2607.19084v1.

All evidence uses rational arithmetic or symbolic polynomial identities. Infinite
sequences are finite sums of geometric sequences; no truncation defines them.
The inequalities certify the whole interval 0 < delta <= 1/1000. Root-limit
and nonconvergence implications are explained in README.md, not numerically
approximated here. Requires Python 3 and SymPy.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import sympy as S


def equal(left, right):
    assert S.cancel(left - right) == 0, (left, right)


def sym(q):
    q = Q(q)
    return S.Rational(q.numerator, q.denominator)


p, eps, r, lam = Q(3, 4), Q(1, 2), Q(4, 5), Q(73, 80)
q = 1 - p
u = [(Q(7, 6), Q(1, 2)), (Q(28, 345), r), (-Q(70, 69), Q(5, 12))]
v = [(-Q(153, 286), Q(1, 2)), (Q(18972, 1495), -r), (-Q(1785, 253), -Q(5, 12))]
z = S.Symbol('z')
delta = S.Symbol('delta', real=True)


def coeff(terms, k):
    return sum((c * s ** k for c, s in terms), Q(0))


def gen(terms):
    return sum((sym(c) / (1 - sym(s) * z) for c, s in terms), S.Integer(0))


def linear_gen(A, x0):
    return sym(p) * (A - x0) / z + sym(q) * z * A + sym(p * (1 - eps)) * x0 / (2 - z)


# Fidelity: reproduce the source's general normalization identity and QSD
# generating-function formula (Section 4.1, label eq:Az), independently of u,v.
P, E, alpha = S.symbols('P E alpha')
L = 1 - P * E * alpha
source_A = P * alpha * (2 - (2 - E) * z) / ((2 - z) * ((1 - P) * z ** 2 - L * z + P))
source_M = P * (source_A - alpha) / z + (1 - P) * z * source_A + P * alpha * (1 - E) / (2 - z)
equal(source_A.subs(z, 0), alpha)
equal(source_A.subs(z, 1), 1)
equal(source_M, L * source_A)
A1, x0 = S.symbols('A1 x0')
equal(P * (A1 - x0) + (1 - P) * A1 + P * x0 * (1 - E), A1 - P * E * x0)
print('FIDELITY_OK: source master equation, mass law, and QSD generating function')

U, V = gen(u), gen(v)
expected_U = S.Rational(7, 30) * (1 - S.Rational(3, 4) * z) / ((1 - z / 2) * (1 - S.Rational(4, 5) * z) * (1 - S.Rational(5, 12) * z))
expected_V = S.Rational(51, 10) * (1 - S.Rational(3, 4) * z) / ((1 - z / 2) * (1 + S.Rational(4, 5) * z) * (1 + S.Rational(5, 12) * z))
equal(U, expected_U)
equal(V, expected_V)
equal(U, source_A.subs({P: sym(p), E: sym(eps), alpha: S.Rational(7, 30)}))
for terms in [u, v]:
    assert all(abs(s) < 1 for _, s in terms)
    assert sum((c / (1 - s) for c, s in terms), Q(0)) == 1
assert coeff(u, 0) == Q(7, 30)
assert coeff(v, 0) == Q(51, 10)
print('GENERATING_FUNCTIONS_OK: U(1)=V(1)=1, u0=7/30, v0=51/10')

# All-index eigenidentities: the coefficient of each independent geometric
# mode in T(x)-eigenvalue*x vanishes; the reset contributes only at s=1/2.
for terms, eigenvalue in [(u, lam), (v, -lam)]:
    for c, s in terms:
        residual = c * (eigenvalue - p * s - q / s)
        target = (p * (1 - eps) / 2) * coeff(terms, 0) if s == Q(1, 2) else Q(0)
        assert residual == target
    assert p * coeff(terms, 1) + p * (1 - eps) / 2 * coeff(terms, 0) == eigenvalue * coeff(terms, 0)
    equal(linear_gen(gen(terms), sym(coeff(terms, 0))), sym(eigenvalue) * gen(terms))
assert coeff(u, 1) == Q(203, 900)
assert coeff(v, 1) == -Q(187, 25)
print('EIGENIDENTITIES_OK: Tu=(73/80)u and Tv=-(73/80)v for every index')

# Uniform bounds for all k >= 0: (1/2)^k >= (5/12)^k and all |s|<=r.
assert Q(1, 2) >= Q(5, 12) > 0
assert Q(7, 6) - Q(70, 69) == Q(7, 46) > 0
lower_u = Q(28, 345)
abs_u = sum((abs(c) for c, _ in u), Q(0))
abs_v = sum((abs(c) for c, _ in v), Q(0))
assert abs_u < 3
assert abs_v == Q(60639, 2990) < 21
assert all(abs(s) <= r for _, s in u + v)
max_delta = Q(1, 1000)
margin = lower_u - 21 * max_delta
assert margin == Q(4151, 69000) > 0
assert (3 + 21 * max_delta) / (1 - max_delta) < 4
print('POSITIVE_FAMILY_OK: all 0<delta<=1/1000, margin=4151/69000, common ratio=4/5')

# Normalized full rational functions, with symbolic delta (the interval above
# guarantees all denominators are positive). These are whole-sequence checks.
a = (U + delta * V) / (1 + delta)
b = (U - delta * V) / (1 - delta)
a0 = a.subs(z, 0)
b0 = b.subs(z, 0)
Da = 1 - sym(p * eps) * a0
Db = 1 - sym(p * eps) * b0
equal(a.subs(z, 1), 1)
equal(b.subs(z, 1), 1)
equal(Da, sym(lam) * (1 - delta) / (1 + delta))
equal(Db, sym(lam) * (1 + delta) / (1 - delta))
equal(linear_gen(a, a0) / Da, b)
equal(linear_gen(b, b0) / Db, a)
equal(a0 - b0, S.Rational(146, 15) * delta / (1 - delta ** 2))
print('TWO_CYCLE_OK: F(a_delta)=b_delta, F(b_delta)=a_delta, a0-b0=146*delta/(15*(1-delta**2))')

sample = sym(max_delta)
for expr, target in [(a0, Q(7153, 30030)), (b0, Q(6847, 29970)),
                     (a0 - b0, Q(29200, 2999997)),
                     (Da, Q(72927, 80080)), (Db, Q(73073, 79920))]:
    equal(expr.subs(delta, sample), sym(target))
# Small-index recurrence checks independently evaluate the defining operator.
for k in range(101):
    def T(terms):
        return (p * coeff(terms, k + 1)
                + (q * coeff(terms, k - 1) if k else Q(0))
                + p * (1 - eps) * coeff(terms, 0) * Q(1, 2) ** (k + 1))
    assert T(u) == lam * coeff(u, k)
    assert T(v) == -lam * coeff(v, k)
    assert coeff(u, k) >= lower_u * r ** k
    assert abs(coeff(v, k)) <= 21 * r ** k
    ak = (coeff(u, k) + max_delta * coeff(v, k)) / (1 + max_delta)
    bk = (coeff(u, k) - max_delta * coeff(v, k)) / (1 - max_delta)
    assert ak > 0 and bk > 0
assert p > 0 and q > 0 and p * (1 - eps) / 2 == Q(3, 16) > 0
print('INDEPENDENT_COORDINATE_CHECK_OK: k=0..100; self-loop=3/16; both adjacent transitions positive')
print('SAMPLE: a0=7153/30030, b0=6847/29970, gap=29200/2999997')
print('SAMPLE_SURVIVAL: 72927/80080, 73073/79920')
print('CHECK_SHA256:', sha256(Path(__file__).read_bytes()).hexdigest())
print('ALL_CHECKS_OK')
