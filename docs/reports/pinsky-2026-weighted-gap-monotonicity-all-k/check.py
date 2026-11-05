#!/usr/bin/env python3
"""Exact checks for Pinsky arXiv:2608.04730v1, equation (1.3).

All bounds use integers or rational arithmetic. No numerical quadrature or floating
point enters the checks. The finite window does not prove the universal claim.

Write H = sum(1/j, 1<=j<k) and exp(2 sum(s**j/j)) = sum(c_n*s**n).
Then c_0=1 and n*c_n=2*sum(c_{n-j}, 1<=j<k, j<=n).
Positive coefficients permit a rigorous tail estimate using their total exp(2H).
Directed integer rounding encloses every c_n in units of SCALE.
"""
from fractions import Fraction as F
from math import factorial
import sympy as S

SCALE = 10**60
TRUNCATION = 2000


def ceildiv(a, b):
    return -(-a // b)


def exp_bounds(x):
    """Taylor lower bound and geometric majorant for the positive remainder."""
    assert x >= 0
    term = F(1)
    lower = term
    for n in range(1, 201):
        term *= x / n
        lower += term
    first = term * x / 201
    ratio = x / 202
    assert ratio < 1
    return lower, lower + first / (1 - ratio)


def density_bounds(k, last_l):
    h = sum((F(1, j) for j in range(1, k)), F(0))
    exp_lo, exp_hi = exp_bounds(2 * h)
    elo = exp_lo.numerator * SCALE // exp_lo.denominator
    ehi = ceildiv(exp_hi.numerator * SCALE, exp_hi.denominator)
    lo, hi = [SCALE], [SCALE]
    window_lo = window_hi = SCALE
    for n in range(1, TRUNCATION + 1):
        lo.append(2 * window_lo // n)
        hi.append(ceildiv(2 * window_hi, n))
        window_lo += lo[n]
        window_hi += hi[n]
        if n >= k - 1:
            window_lo -= lo[n - k + 1]
            window_hi -= hi[n - k + 1]
    tail = max(0, ehi - sum(lo))
    out = []
    for l in range(last_l + 1):
        lower = sum(2 * a // ((n+l+1)*(n+l+2)) for n, a in enumerate(lo))
        upper = sum(ceildiv(2*a, (n+l+1)*(n+l+2)) for n, a in enumerate(hi))
        upper += ceildiv(2*tail, (TRUNCATION+l+2)*(TRUNCATION+l+3))
        # Integral numerator and exp(2H) have the same scale, which cancels.
        out.append((F(lower, ehi), F(upper, elo)))
    return out


def factorial_moment(poly, t):
    return S.expand(sum(c * S.factorial(m[0]) for m, c in S.Poly(poly, t).terms()))


def polynomial_certificate():
    t, u = S.symbols('t u', real=True)
    p = (6 - 18*t + 9*t**2 - t**3) / 34
    q = (24 - 36*t + 12*t**2 - t**3) / 73
    low, r1 = S.div(1 - p**2, 1+t, t)
    high, r2 = S.div(1 + t*q**2, 1+t, t)
    assert r1 == r2 == 0
    assert S.expand(1 - (1+t)*low - p**2) == 0
    assert S.expand((1+t)*high - 1 - t*q**2) == 0
    assert factorial_moment(low, t) == S.Rational(10, 17)
    assert factorial_moment(high, t) == S.Rational(44, 73)
    assert sum((F(1, factorial(j)) for j in range(5)), F(0)) == F(65, 24)
    c_margin = S.Rational(229, 288) - S.Rational(1056, 4745) - S.Rational(57, 100)
    assert c_margin > 0
    z = S.Rational(57, 50) - 2*u + u**2/2 - u**3/9
    head = S.integrate(u**2 * sum(z**j/S.factorial(j) for j in range(4)), (u, 0, 1))
    assert head == S.Rational(39235476059, 126299250000)
    bound = 1 - head + 2*S.Rational(24, 65)*(1 - S.Rational(10, 17))
    nominated = S.Rational(414358946399, 416821204800)
    assert bound == S.Rational(27728403640961, 27912134250000)
    assert bound < nominated < 1
    print('E1_POLYNOMIAL_LOWER_MOMENT=10/17')
    print('E1_POLYNOMIAL_UPPER_MOMENT=44/73')
    print('HEAD_INTEGRAL_LOWER=' + str(head))
    print('INDEPENDENT_SCALAR_BOUND=' + str(bound))
    print('NOMINATED_SCALAR_BOUND=' + str(nominated))
    print('NOMINATED_BOUND_MARGIN=' + str(1-nominated))


def fidelity():
    s = S.symbols('s', real=True)
    k2 = S.integrate(2*(1-s)*s*S.exp(2*(s-1)), (s, 0, 1))
    assert S.simplify(k2 - S.exp(-2)) == 0
    print('SOURCE_FIDELITY_g_2_1=exp(-2)')
    for k in range(2, 61):
        # Divide F' = -2(1-s)s^(k-1)exp(2 sum s^j/j) by exp(...).
        pprime = 2*sum(s**j for j in range(k-1))
        residual = -2*(1-s) + (1-s)**2*pprime + 2*(1-s)*s**(k-1)
        assert S.Poly(S.expand(residual), s).is_zero
    print('SOURCE_LONGEST_GAP_DERIVATIVE_IDENTITIES=59')


def main():
    fidelity()
    polynomial_certificate()
    minimum = None
    stronger_minimum = None
    checked = 0
    for k in range(2, 61):
        bounds = density_bounds(k, k-1)
        for l in range(k-1):
            gl, gu = bounds[l]
            nl, nu = bounds[l+1]
            dlo, dhi = (l+1)*nl-l*gu, (l+1)*nu-l*gl
            strong = dlo-gu/F(2*k-1)
            assert dlo > 0, (k, l, dlo)
            assert strong > 0, (k, l, strong)
            checked += 1
            if minimum is None or dlo < minimum[0]:
                minimum = (dlo, dhi, k, l)
            if stronger_minimum is None or strong < stronger_minimum[0]:
                stronger_minimum = (strong, k, l)
        if k == 10:
            gl, gu = bounds[5]
            assert F(5588, 1000000) < gl < gu < F(5590, 1000000)
            print('TABLE1_CORRECTION_g_10_5_IN=(5588/1000000,5590/1000000)')
    dlo, dhi, k, l = minimum
    # Short decimal *rational* enclosure; no float conversion.
    report_scale = 10**18
    lower = dlo.numerator*report_scale//dlo.denominator
    upper = ceildiv(dhi.numerator*report_scale, dhi.denominator)
    assert lower > 0
    # Certify the location as well as a positive lower bound.
    for kk in range(2, 61):
        if kk == k:
            continue
        bb = density_bounds(kk, kk-1)
        assert all((ll+1)*bb[ll+1][0]-ll*bb[ll][1] > dhi for ll in range(kk-1))
    bb = density_bounds(k, k-1)
    assert all((ll+1)*bb[ll+1][0]-ll*bb[ll][1] > dhi for ll in range(k-1) if ll != l)
    print('WINDOW=2<=k<=60,0<=l<=k-2')
    print('CHECKED_INEQUALITIES=' + str(checked))
    print('MINIMUM_SLACK_LOCATION=k:' + str(k) + ',l:' + str(l))
    print('MINIMUM_SLACK_RATIONAL_ENCLOSURE=[' + str(lower) + '/' + str(report_scale) + ',' + str(upper) + '/' + str(report_scale) + ']')
    stronger, sk, sl = stronger_minimum
    slb = stronger.numerator*report_scale//stronger.denominator
    print('STRONGER_MARGIN_MINIMUM_LOWER_BOUND=' + str(slb) + '/' + str(report_scale))
    print('STRONGER_MARGIN_MINIMUM_LOWER_BOUND_LOCATION=k:' + str(sk) + ',l:' + str(sl))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
