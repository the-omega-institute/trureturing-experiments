#!/usr/bin/env python3
"""Exact checks for Karp's SF21 Conjecture 4; no floating-point evidence.

Direct rational series enclosures are independent of the endpoint expansion.
The symbolic checks verify the algebra used in the proposed asymptotic family;
they do not certify the analytic continuation formula or its uniform remainder.
"""
from fractions import Fraction as Q
import hashlib
from pathlib import Path
import sympy as sp


def ceil_div(a, b):
    assert b > 0
    return -((-a) // b)


def series_enclosure(nu, y, nterms, bits=192):
    """Enclose sum c_j*y**j using outward-rounded dyadic intervals.

    For j >= nterms > nu > 0, the exact term ratio is
      y*(1 - 1/(j+1) - nu*(nu+1)/(j+1)**2),
    which lies in (0,y). Thus |tail| <= |term_nterms|*y/(1-y).
    Integer floor/ceiling at each operation preserves enclosure, even for
    negative early coefficients. Rounding is exact integer arithmetic.
    """
    assert isinstance(nu, Q) and isinstance(y, Q)
    assert 0 < nu < nterms and 0 < y < 1
    scale = 1 << bits
    low = high = scale
    sum_low = sum_high = scale
    p, q = nu.numerator, nu.denominator
    r, s = y.numerator, y.denominator
    for j in range(nterms):
        numerator = (j*q - p) * ((j+1)*q + p) * r
        denominator = q*q * (j+1)*(j+1) * s
        a, b = low*numerator, high*numerator
        low, high = min(a,b)//denominator, ceil_div(max(a,b),denominator)
        sum_low += low
        sum_high += high
    ratio0 = (Q(nterms)-nu)*(Q(nterms)+nu+1)/(nterms+1)**2
    assert 0 < ratio0 < 1
    tail = Q(max(abs(low),abs(high)),scale) * y/(1-y)
    return Q(sum_low,scale)-tail, Q(sum_high,scale)+tail


def mul_interval(a,b):
    p = [x*y for x in a for y in b]
    return min(p), max(p)


def square_interval(a):
    lo,hi = a
    return (Q(0) if lo <= 0 <= hi else min(lo*lo,hi*hi), max(lo*lo,hi*hi))


def coarse(a, scale=10**8):
    lo,hi = a
    return Q((lo*scale).__floor__(),scale), Q((hi*scale).__ceil__(),scale)


def f_poly(n,y):
    return sp.expand(sum(sp.rf(-n,j)*sp.rf(n+1,j)*y**j/sp.factorial(j)**2
                         for j in range(n+1)))


def fidelity():
    y = sp.Symbol('y')
    assert sp.expand(f_poly(1,y)**2-f_poly(0,y)*f_poly(2,y)-2*y*(1-y)) == 0
    print('FIDELITY: source integer case mu=1, delta=1: T=2*y*(1-y)>0 for 0<y<1')
    count = 0
    for n in range(1,13):
        expr = f_poly(n,y)**2 - f_poly(n-1,y)*f_poly(n+1,y)
        for denominator in range(2,14):
            for numerator in range(1,denominator):
                assert expr.subs(y,sp.Rational(numerator,denominator)) > 0
                count += 1
    print(f'FIDELITY: exact positive integer-degree grid checks={count}')
    # Check the interval implementation against exact terminating series.
    for n in range(1,7):
        for y0 in [Q(1,4), Q(1,2), Q(3,4)]:
            lo,hi = series_enclosure(Q(n),y0,16)
            exact = Q(f_poly(n,sp.Rational(y0.numerator,y0.denominator)))
            assert lo <= exact <= hi
    print('ENCLOSURE CROSS-CHECK: 18 terminating-series exact values enclosed')


def family_algebra():
    m,b,a,t = sp.symbols('m b a t')
    bm = b+2/m
    bp = b-2/(m+1)
    assert sp.factor(2*b-bm-bp+2/(m*(m+1))) == 0
    leading = sp.expand(a*a*(t+b)**2-a*a*(t+bm)*(t+bp))
    assert sp.factor(leading.coeff(t,1)+2*a*a/(m*(m+1))) == 0
    assert leading.coeff(t,2) == 0
    print('FAMILY ALGEBRA: neighboring-digamma difference = -2/(mu*(mu+1))')
    print('FAMILY ALGEBRA: delta=1 linear coefficient = -2*a_mu**2/(mu*(mu+1))')
    d = sp.Symbol('d', real=True)
    mreal = sp.Symbol('mreal', real=True)
    trig = sp.expand_trig(sp.sin(sp.pi*(mreal-d))*sp.sin(sp.pi*(mreal+d)))
    assert sp.trigsimp(trig-sp.sin(sp.pi*mreal)**2+sp.sin(sp.pi*d)**2) == 0
    assert sp.simplify(2*sp.sin(sp.pi*sp.Rational(3,2))**2 /
                       (sp.pi**2*sp.Rational(3,2)*sp.Rational(5,2))) == 8/(15*sp.pi**2)
    z = sp.Symbol('z')
    assert sp.series(sp.sin(sp.pi*z)**2/sp.pi**2,z,0,6).removeO().coeff(z,2) == 1
    print('FAMILY ALGEBRA: C_(3/2)=8/(15*pi**2)>0; quadratic moving-shift coefficient=1')
    print('ANALYTIC BOUNDARY: DLMF 15.8.10 and the uniform endpoint remainder are prose obligations')


def finite_counterexamples():
    mu = Q(3,2)
    # Rational points near the joint endpoint; no gamma/digamma evaluation.
    witnesses = [(Q(127,128), 8), (Q(63,64), 10), (Q(63,64), 12)]
    for delta,k in witnesses:
        y = 1-Q(1,2**k)
        nterms = 32*(2**k)
        assert 0 < delta < 1 and 0 < y < 1 and mu-delta > 0
        center = series_enclosure(mu,y,nterms)
        minus = series_enclosure(mu-delta,y,nterms)
        plus = series_enclosure(mu+delta,y,nterms)
        center_sq = square_interval(center)
        product = mul_interval(minus,plus)
        turan = center_sq[0]-product[1], center_sq[1]-product[0]
        lo,hi = coarse(turan)
        print(f'WITNESS mu={mu} delta={delta} y={y} nterms={nterms}: {lo} <= T <= {hi}')
        print(f'  F_minus in {coarse(minus)}; F_center in {coarse(center)}; F_plus in {coarse(plus)}')
        assert hi < 0, 'proposed strict counterexample did not certify'
    print(f'REFUTATION CHECK: {len(witnesses)} exact rational strict counterexamples; every degree is positive')


if __name__ == '__main__':
    print('CHECK_SHA256='+hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ARITHMETIC: integers, fractions.Fraction, sympy exact algebra; no evalf or float')
    fidelity()
    family_algebra()
    finite_counterexamples()
    print('ALL_CHECKS_OK')
