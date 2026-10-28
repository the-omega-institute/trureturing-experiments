#!/usr/bin/env python3
"""Exact source-fidelity and interval certificates for YCZ Conjecture 1(ii).

Analytic inputs (not proved by running this program): DLMF 5.7.6,
termwise derivatives on compact subsets of u > 0, gamma < 1, the
atanh series for log, and the decreasing-function integral tail bound.
All numerical comparisons and symbolic identities below are exact.
"""
from fractions import Fraction as Q
import sympy as s


def main():
    x, a, n, m = s.symbols('x a n m', real=True)
    A, eps, r = Q(155, 294), Q(1, 10000), Q(1, 64)
    u = r + Q(1, 2)
    assert 0 < A - eps < A and 0 < r and 0 < u
    assert Q(7, 40) * A == Q(31, 336)
    print('DOMAIN: 0 < a_low =', A - eps, '<= a <=', A, '; x =', r)

    # Source: Remark 4, rational improvement identity following Corollary 4.
    z = n + s.Rational(1, 2)
    new = (z**2 + s.Rational(2071, 5880)) / (24*z**2*(z**2+s.Rational(155,294)))
    old = (z**2-s.Rational(7,40))/(24*(z**4-s.Rational(31,336)))
    source_rhs = -s.Rational(64201,120) / ((2*n+1)**2*(588*n**2+588*n+457)*
                                        (168*n**4+336*n**3+252*n**2+84*n-5))
    assert s.cancel(new-old-source_rhs) == 0
    # For n >= 1 every denominator factor is positive: shift n = m+1.
    for poly in [2*n+1,588*n**2+588*n+457,168*n**4+336*n**3+252*n**2+84*n-5]:
        assert all(c > 0 for c in s.Poly(s.expand(poly.subs(n,m+1)),m).all_coeffs())
    assert all(source_rhs.subs(n,k) < 0 for k in range(1,33))
    print('FIDELITY: source Remark 4 identity verified symbolically; sign for every integer n >= 1; n=1..32 cross-check')

    # Independent R lower bound: psi(u)=-gamma-1/u+sum_{k>=1}(1/k-1/(k+u)).
    # gamma < 1.  log 2 > 2*(1/3 + (1/3)^3/3) = 56/81.
    # Discard the positive tail after k=3; -log(r)=6 log 2.
    log2_lower = 2*(Q(1,3)+Q(1,3)**3/3)
    psi_lower = -1-1/u+sum((1/Q(k)-1/(Q(k)+u) for k in range(1,4)), Q(0))
    R_lower = psi_lower+6*log2_lower
    assert log2_lower == Q(56,81)
    assert R_lower == Q(394328567,231912450) and R_lower > Q(5,3)
    print('R(r) >', R_lower, '> 5/3')

    # psi'(u)=sum_{k>=0}(k+u)^(-2)>0, so R'(r)>-1/r.
    # -psi''(u)=2 sum_{k>=0}(k+u)^(-3).
    # k>=2 tail < integral_1^infty (u+t)^(-3)dt = 1/(2(u+1)^2).
    psi2_abs_upper = 2*(1/u**3+1/(u+1)**3+1/(2*(u+1)**2))
    assert psi2_abs_upper == Q(511623041024,32798729601)
    assert psi2_abs_upper < 16 < 18
    print("R'(r) >", -1/r)
    print("-psi''(u) <", psi2_abs_upper, "< 16 < 18; R''(r) >",1/r**2-18)

    c = s.Rational(7,40)*(a-s.Rational(A.numerator,A.denominator))
    P = x**4+a*x**2+c
    RR = s.Function('R')
    G = 24*P*RR(x)-(x**2-s.Rational(7,40)+a)
    formula = 24*((12*x**2+2*a)*RR(x)+(8*x**3+4*a*x)*s.diff(RR(x),x)+P*s.diff(RR(x),x,2))-2
    assert s.expand(s.diff(G,x,2)-formula) == 0
    reduced = 24*(13*x**2+a/3+c/x**2-18*P)-2
    substituted = formula.subs({s.diff(RR(x),x,2):1/x**2-18,
                               s.diff(RR(x),x):-1/x,RR(x):s.Rational(5,3)})
    assert s.cancel(substituted-reduced) == 0
    print("DERIVATIVE: G_a'' product rule and substitution identities verified symbolically")

    Pmin = r**4+(A-eps)*r**2-Q(7,40)*eps
    Pmax = r**4+A*r**2
    assert Pmin == Q(857398127,7707033600000) and Pmin > 0
    assert Pmax == Q(317587,2466250752) and Pmin <= Pmax
    assert 12*r**2+2*(A-eps) > 0 and 8*r**3+4*(A-eps)*r > 0
    print('INTERVAL: P(r) in [',Pmin,',',Pmax,']; all three multiplier coefficients positive')

    # Uniform inequalities use a>=A-eps, c>=-7eps/40 and P(r)<=Pmax.
    bound = 24*(13*r**2+(A-eps)/3-Q(7,40)*eps/r**2-18*Pmax)-2
    assert bound == Q(249085992523,481689600000) and bound > Q(1,2)
    # Independent affine endpoint check, keeping dependence c=7(a-A)/40 intact.
    sr = s.Rational(r.numerator,r.denominator)
    aa = s.Rational((A-eps).numerator,(A-eps).denominator)
    exact_affine = s.cancel(reduced.subs(x,sr))
    slope = s.diff(exact_affine,a)
    assert slope > 0 and s.diff(exact_affine,a,2) == 0
    assert exact_affine.subs(a,aa) >= s.Rational(bound.numerator,bound.denominator)
    print("UNIFORM: G_a''(1/64) >",bound,"> 1/2 for every a in [155/294-1/10000,155/294]")
    print("SETTLEMENT CHECK: (-G_a)''(1/64) < -1/2; violates the n=2 CM condition in the conjectured sufficient direction")
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
