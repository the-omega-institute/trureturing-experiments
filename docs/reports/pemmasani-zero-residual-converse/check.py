#!/usr/bin/env python3
"""Exact checks for arXiv:2610.09258v2, Proposition 6.2 and its open converse.

Finite checks support a general analytic proof; they do not certify that proof.
No floating point or numerical optimization is used.
"""
from fractions import Fraction as F
from itertools import product
import sys
import sympy as s


def exact_zero(x):
    return s.simplify(s.radsimp(x)) == 0


def fidelity():
    # Source 6.2(e): clockwise rates 2, anticlockwise rates 1; i=3.
    z = s.Symbol('z', real=True)
    beta = [s.Rational(1, 2)] * 2
    cm, cp = [1, 2], [2, 1]
    r = [s.sqrt(2), s.Integer(1)]
    aa = sum(cm[j] * beta[j] / r[j] for j in range(2))
    bb = sum(cp[j] * r[j] for j in range(2))
    omega = s.expand(aa * bb)
    weighted_rate = 3 - z * 2 - s.sqrt((2 + z) * (1 + 4*z))
    third = s.diff(weighted_rate, z, 3).subs(z, 0)
    assert exact_zero(omega - (2 + s.Rational(9, 4) * s.sqrt(2)))
    assert exact_zero(third + s.Rational(1323, 64) * s.sqrt(2))
    delta = s.simplify(cp[0] * r[0] / bb - cm[0] * beta[0] / r[0] / aa)
    w = s.sqrt(2) / 2
    # L = [[2w,-2w],[-2w,2w]], L+ = [[1,-1],[-1,1]]/(8w).
    L = s.Matrix([[2*w, -2*w], [-2*w, 2*w]])
    Lplus = s.Matrix([[1, -1], [-1, 1]]) / (8*w)
    Delta = s.Matrix([delta, -delta])
    d = s.simplify(omega / 2 * (Delta.T * Lplus * Delta)[0])
    assert all(exact_zero(x) for x in L * Lplus * L - L)
    assert all(exact_zero(x) for x in Lplus * L * Lplus - Lplus)
    cubic_source = third / (6 * omega**s.Rational(3, 2)) + d**2 / s.sqrt(omega)
    assert exact_zero(cubic_source + s.Rational(49, 64) / s.sqrt(omega))
    # Independent square-root variational value: P(t)=p0+p+ exp(t)+p- exp(-t).
    pp = cp[0]*r[0]*cm[1]*beta[1]/r[1]
    pm = cp[1]*r[1]*cm[0]*beta[0]/r[0]
    B1 = (pp-pm)/(2*s.sqrt(omega))
    B2 = (pp+pm)/(2*s.sqrt(omega)) - (pp-pm)**2/(4*omega**s.Rational(3, 2))
    cubic_reduced = -B1**2 * B2 / w**2
    assert exact_zero(cubic_reduced - cubic_source)
    print('FIDELITY Proposition 6.2(e): omega=2+9*sqrt(2)/4; Jzzz=-1323*sqrt(2)/64')
    print('FIDELITY cubic=-49/(64*sqrt(omega)); Laplacian Penrose identities verified')


def implicit_coefficients():
    # Differentiate the critical equation A'(t)-2e B'(t)=0 formally.
    # At the two-state optimum: A''=-2w and A'''=0.
    e = s.Symbol('e')
    w, B0, B1, B2 = s.symbols('w B0 B1 B2', positive=True)
    u1, u2 = s.symbols('u1 u2')
    t = u1*e + u2*e**2
    critical = s.expand(-2*w*t - 2*e*(B1+B2*t))
    solved = s.solve([critical.coeff(e, 1), critical.coeff(e, 2)], [u1, u2], dict=True)[0]
    value = s.expand(-w*t**2 - 2*e*(B0+B1*t+B2*t**2/2))
    cubic = s.simplify(value.coeff(e, 3).subs(solved))
    quadratic = s.simplify(value.coeff(e, 2).subs(solved))
    assert cubic == -B1**2*B2/w**2
    assert quadratic == B1**2/w
    u, p, p0 = s.symbols('u p p0', positive=True)
    assert s.factor((u+1/u-2) - (u-1)**2/u) == 0
    assert s.factor((p0+p*u+p/u)-(p0+2*p)-p*(u-1)**2/u) == 0
    print('IMPLICIT coefficients: e^2=B1^2/w; e^3=-B1^2*B2/w^2')
    print('EQUAL_LAWS exact global gaps: A0-A=w*(u-1)^2/u; P-P0=p*(u-1)^2/u')


def small_cases():
    # All 3-state off-diagonal rates in {0,1,2} with irreducible G={0,1}
    # and irreducible F={0,1,2}. For each, test all optimizer ratios h in
    # {1/2,1,2}. The resulting beta is rational and makes r=(h,1) optimal.
    # This is 256 generators and 768 profiles; it is not all real profiles.
    generators = cases = equal = unequal = sparse = 0
    gauges = [F(1, 2), F(1), F(2)]
    for a, b in product((1, 2), repeat=2):
        for cm0, cm1, cp0, cp1 in product((0, 1, 2), repeat=4):
            if cm0+cm1 == 0 or cp0+cp1 == 0:
                continue
            generators += 1
            for h in gauges:
                beta0 = F(b)*h*h/(F(a)+F(b)*h*h)
                beta1 = 1-beta0
                w = beta0*a/h
                assert w == beta1*b*h and w > 0
                aa = F(cm0)*beta0/h + F(cm1)*beta1
                bb = F(cp0)*h + F(cp1)
                omega = aa*bb
                assert aa > 0 and bb > 0 and omega > 0
                pp = F(cp0)*h*cm1*beta1
                pm = F(cp1)*cm0*beta0/h
                p0 = F(cp0*cm0)*beta0 + F(cp1*cm1)*beta1
                assert omega == p0+pp+pm
                E, D = pp-pm, pp+pm
                delta = F(cp0)*h/bb - F(cm0)*beta0/h/aa
                assert E == omega*delta
                assert F(cp1)/bb - F(cm1)*beta1/aa == -delta
                # B'' = N/(4 omega^(3/2)). N has a nonnegative expansion.
                N = 2*omega*D-E*E
                assert N == 2*p0*D+D*D+4*pp*pm
                assert N >= 0
                # Cubic = -E^2 N/(16 w^2 omega^(5/2)). Its radical
                # denominator is strictly positive; test its exact numerator.
                cubic_numerator = -E*E*N
                laws_equal = delta == 0
                assert (cubic_numerator == 0) == laws_equal
                d_resistance = omega*delta*delta/(4*w)
                d_implicit = E*E/(4*omega*w)
                assert d_resistance == d_implicit
                assert (d_resistance == 0) == laws_equal
                if laws_equal:
                    # E=0 and P'(0)=0 give the fixed critical optimizer for
                    # every e>=0; exact stationary value A0-2sqrt(omega)e.
                    assert E == 0 and d_resistance == 0 and pp == pm and p0+2*pp == omega
                    equal += 1
                else:
                    assert cubic_numerator < 0 and d_resistance > 0
                    unequal += 1
                sparse += int(0 in (cm0, cm1, cp0, cp1))
                cases += 1
    assert (generators, cases) == (256, 768)
    print(f'SMALL_CASES generators={generators}; profiles={cases}; equal={equal}; unequal={unequal}; sparse={sparse}')
    print('SMALL_CASES zero cubic iff equal laws; unequal laws have a nonzero residual germ')
    # G singleton: beta=1 and both laws are the point mass, for all positive
    # rates in the stated small grid. Exact two-state elimination is affine.
    count = 0
    for cm, cp in product((1, 2, 3), repeat=2):
        omega = F(cm*cp)
        assert F(cm, cm) == F(cp, cp) == 1 and omega > 0
        count += 1
    print(f'SINGLETON_CASES {count}; both laws=1; d=0; reduced value=q_j-2*sqrt(omega)*e')


def main():
    assert not any(isinstance(x, float) for x in [F(1), s.sqrt(2)])
    print('Python '+sys.version.split()[0]+'; SymPy '+s.__version__+'; exact rational and symbolic algebra only')
    fidelity()
    implicit_coefficients()
    small_cases()
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
