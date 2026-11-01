#!/usr/bin/env python3
"""Exact checks for KVZ arXiv:2408.01755v2, Conjecture 1, first sentence.

The finite grid is a diagnostic, not a proof of the universal assertion.
Every interval encloses the infinite positive series, including its tail.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import sympy as s


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


ORDERS = tuple(map(Q, ['-3/4', '-1/2', '-1/4', '0', '1/2', '1', '3/2', '2']))
SCALES = ((Q(1, 2), Q(1)), (Q(1), Q(1)), (Q(1), Q(2)), (Q(2), Q(2)))
POINTS = tuple(map(Q, ['1/8', '1/4', '1/2', '1', '2', '4', '8', '16']))
N = 128


def coefficients(nu, a, count):
    """I_nu(a*x) / [(a/2)^nu/Gamma(nu+1)*x^nu], in powers x^(2k)."""
    out = [Q(1)]
    for k in range(count - 1):
        out.append(out[-1] * a*a / (4 * (k+1) * (nu+1+k)))
    return out


# Fidelity: Section 4's equal-order quotient is decreasing (weakly if scales equal).
# Its normalized coefficient ratios are exactly (a1/a2)^(2k).
fidelity_coefficients = 0
for nu, (a1, a2) in product(ORDERS, SCALES):
    p = coefficients(nu, a1, 18)
    q = coefficients(nu, a2, 18)
    ratios = [u/v for u, v in zip(p, q)]
    require(ratios == [(a1/a2)**(2*k) for k in range(18)], 'coefficient-ratio fidelity')
    require(all(v <= u for u, v in zip(ratios, ratios[1:])), 'equal-order coefficient signs')
    fidelity_coefficients += len(ratios)
print('FIDELITY: equal-order coefficient ratios, exact identities =', fidelity_coefficients)

# ODE recurrence and contiguous identities, all rational cases through degree 17.
recurrence_checks = 0
for nu, a in product(ORDERS, (Q(1, 2), Q(1), Q(2))):
    c = coefficients(nu, a, 18)
    for k in range(1, 18):
        require(((nu+2*k)**2-nu**2)*c[k] == a*a*c[k-1], 'modified-Bessel ODE recurrence')
        recurrence_checks += 1
for nu in ORDERS:
    b = nu+1
    c0 = coefficients(nu, Q(2), 18)  # F_b(t), since a^2/4 = 1
    c1 = coefficients(nu+1, Q(2), 18)
    c2 = coefficients(nu+2, Q(2), 18)
    for k in range(17):
        require((k+1)*c0[k+1] == c1[k]/b, 'first derivative coefficient')
        require(c0[k] == c1[k] + (c2[k-1]/(b*(b+1)) if k else 0), 'contiguous coefficient')
        recurrence_checks += 2
print('EXACT SERIES: ODE/derivative/contiguous coefficient checks =', recurrence_checks)

# Symbolic identities, with no numerical evaluations or special-function black box.
x, a1, a2, nu1, nu2 = s.symbols('x a1 a2 nu1 nu2', nonzero=True)
y, z = s.Function('y')(x), s.Function('z')(x)
W = x*(s.diff(y, x)*z-y*s.diff(z, x))
y2 = (a1**2+nu1**2/x**2)*y-s.diff(y, x)/x
z2 = (a2**2+nu2**2/x**2)*z-s.diff(z, x)/x
residual = s.diff(W, x).subs({s.diff(y, x, 2): y2, s.diff(z, x, 2): z2})
target = (nu1**2-nu2**2-(a2**2-a1**2)*x**2)*y*z/x
require(s.simplify(residual-target) == 0, 'weighted-Wronskian identity')
require(s.simplify(s.diff(y/z, x)*x*z**2-W) == 0, 'quotient derivative identity')
P = 1+a1**2*x**2/(4*(nu1+1))
Qpoly = 1+a2**2*x**2/(4*(nu2+1))
H = s.expand((nu1-nu2)*P*Qpoly+x*(s.diff(P, x)*Qpoly-P*s.diff(Qpoly, x)))
require(s.simplify(H.coeff(x, 0)-(nu1-nu2)) == 0, 'unequal-order endpoint coefficient')
require(s.simplify(H.subs(nu1, nu2).coeff(x, 2)-(a1**2-a2**2)/(2*(nu2+1))) == 0,
        'equal-order endpoint coefficient')
print('SYMBOLIC: weighted Wronskian, quotient derivative, both endpoint coefficients = 0 residual')


@lru_cache(maxsize=None)
def F_interval(b, t):
    """L <= sum t^k/(k! (b)_k) <= U, using an exact geometric tail."""
    require(b > 0 and t > 0, 'positive-series domain')
    term = Q(1)
    lower = Q(0)
    for k in range(N):
        lower += term
        term *= t / ((k+1)*(b+k))
    # term is now term_N. Ratios from term_N onward decrease with k.
    rho = t / ((N+1)*(b+N))
    require(0 < rho < 1, 'geometric-tail ratio')
    return lower, lower+term/(1-rho)


@lru_cache(maxsize=None)
def logarithmic_euler_interval(nu, a, xval):
    """Interval for x*y'(x)/y(x), y(x)=I_nu(a*x)."""
    b, t = nu+1, a*a*xval*xval/4
    L, U = F_interval(b, t)
    Lnext, Unext = F_interval(b+1, t)
    multiplier = 2*t/b
    return nu+multiplier*Lnext/U, nu+multiplier*Unext/L


families = evaluated = equal_order_fidelity = constants = 0
patterns = {}
for nu1, nu2 in product(ORDERS, repeat=2):
    if nu1 < nu2:
        continue
    for a1, a2 in SCALES:
        signs = []
        for xval in POINTS:
            if nu1 == nu2 and a1 == a2:
                sign = 0  # Literal identical numerator and denominator, both positive.
                constants += 1
            else:
                L1, U1 = logarithmic_euler_interval(nu1, a1, xval)
                L2, U2 = logarithmic_euler_interval(nu2, a2, xval)
                lo, hi = L1-U2, U1-L2
                require(lo > 0 or hi < 0, 'interval does not decide derivative sign')
                sign = 1 if lo > 0 else -1
            if nu1 == nu2:
                require(sign <= 0, 'source equal-order monotonicity fidelity')
                equal_order_fidelity += 1
            signs.append(sign)
            evaluated += 1
        # The proof route predicts only a positive-to-negative derivative transition.
        nonzero = [v for v in signs if v]
        require(all(not (u < 0 < v) for u, v in zip(nonzero, nonzero[1:])),
                'forbidden negative-to-positive derivative transition')
        key = ''.join('+' if v > 0 else '-' if v < 0 else '0' for v in signs)
        patterns[key] = patterns.get(key, 0)+1
        families += 1
print('FIDELITY: equal-order infinite-series derivative signs =', equal_order_fidelity)
print('GRID: orders =', ','.join(map(str, ORDERS)))
print('GRID: scale pairs =', ';'.join(str(a)+','+str(b) for a, b in SCALES))
print('GRID: positive points =', ','.join(map(str, POINTS)))
print('GRID: parameter families =', families, '; exact derivative signs =', evaluated,
      '; literal constant cases =', constants, '; tail cutoff =', N)
print('GRID: sign patterns =', '; '.join(k+':'+str(v) for k, v in sorted(patterns.items())))
print('SCOPE: exact finite checks and symbolic identities; no universal or kernel proof claimed')
print('ALL_CHECKS_OK')
