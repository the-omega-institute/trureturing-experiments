#!/usr/bin/env python3
"""Exact checks for Haraux 2611.04738v1, Remark 2.2.

No floating point arithmetic is used. Finite checks test the proposed route;
they do not certify the universal sharp-threshold theorem or its endpoint.
"""
from fractions import Fraction as Q
from math import comb, factorial, ceil
import sympy as sp


def moment(k, m):
    """Integral from 0 to 1 of (1-r*r)^k r^(2*m), by expansion."""
    return sum((Q((-1)**j * comb(k, j), 2*(m+j)+1)
                for j in range(k+1)), Q(0))


def j0_partial(z, order):
    return sum(((-1)**k * (z*z/4)**k / factorial(k)**2
                for k in range(order+1)), Q(0))


def cosine_lower(lam, t, order=24):
    """Lower bound for cos(sqrt(lam)*t), from its rational Taylor series.

    The remainder is bounded by the next term once the subsequent absolute
    terms decrease to zero. We check that ratio directly at the cutoff.
    """
    y = lam*t*t
    term, total = Q(1), Q(1)
    for k in range(1, order+1):
        term *= -y / ((2*k-1)*(2*k))
        total += term
    error = abs(term*y/((2*order+1)*(2*order+2)))
    assert y < (2*order+3)*(2*order+4)
    return total-error


def arctan_bounds(x, order=30):
    assert 0 < x < 1 and order % 2 == 0
    s = sum((Q((-1)**k, 2*k+1)*x**(2*k+1)
             for k in range(order+1)), Q(0))
    return s-x**(2*order+3)/(2*order+3), s


# Source fidelity: Remark 2.1, c=3/4 and l=pi.
t, x = sp.symbols('t x', real=True)
u = sp.sin(t/2)*sp.sin(x)
assert sp.simplify(sp.diff(u,t,2)-sp.diff(u,x,2)-sp.Rational(3,4)*u) == 0
assert sp.simplify(u.subs(x,0)) == sp.simplify(u.subs(x,sp.pi)) == 0
from sympy.calculus.util import function_range
assert function_range(sp.sin(t/2), t, sp.Interval(0,2*sp.pi)) == sp.Interval(0,1)
print('FIDELITY: source c=3/4 solution; PDE, boundaries, exact time range [0,1] OK')

# Transmutation coefficients and the supplied beta-integral formula.
for k in range(21):
    for m in range(21):
        rhs = Q(2**(2*k+1)*factorial(k)*factorial(2*m)*factorial(m+k+1),
                factorial(m)*factorial(2*m+2*k+2))
        assert moment(k,m) == rhs
for N in range(31):
    for i in range(1,N+1):
        m = N-i
        coefficient = moment(i-1,m)/(
            2**(2*i-1)*factorial(i-1)*factorial(i)*factorial(2*m))
        assert coefficient == Q(comb(N,i),factorial(2*N))
print('TRANSMUTATION: 441 beta moments and all coefficients through N=30 OK')
q, z = sp.symbols('q z')
for N in range(19):
    p = sum(sp.binomial(N,i)*q**i*(-z*z)**(N-i) for i in range(N+1))
    assert sp.expand(p-(q-z*z)**N) == 0
# The triangular matrices for T_q T_r=T_(q+r), coefficientwise.
for m in range(13):
    for h in range(13-m):
        for i in range(h+1):
            assert comb(m+i,i)*comb(m+h,h-i) == comb(m+h,h)*comb(h,i)
print('BINOMIAL: (q-z^2)^N through N=18; semigroup through degree 24 OK')
for k in range(31):
    assert moment(k,0)/Q(4**k*factorial(k)**2) == Q(1,factorial(2*k+1))
print('WEIGHT: integral J0(sqrt(c)*sqrt(L^2-t^2)) has sine coefficients through k=30 OK')

# Alternating J0 series: terms decrease from k=1 for z<=3.
for z0 in [Q(3),Q(241,100)]:
    assert z0*z0/4 < 4
    assert j0_partial(z0,8) < 0
print('BESSEL: J0(3) <', j0_partial(Q(3),8), '< 0 (alternating upper bound)')
# P_5 is an alternating lower bound on 0<=z<=12/5. Its Bernstein
# coefficients on y=z^2/4 in [0,36/25] are all strictly positive.
power = [Q((-1)**k,factorial(k)**2) for k in range(6)]
bernstein = [sum((power[k]*Q(36,25)**k*Q(comb(i,k),comb(5,k))
                  for k in range(i+1)),Q(0)) for i in range(6)]
assert all(b>0 for b in bernstein)
print('BESSEL: J0 positive on [0,12/5]; J0(241/100)<0; first zero lies in (12/5,241/100)')
print('BESSEL Bernstein coefficients:', ', '.join(map(str,bernstein)))

# Machin's exact identity: pi=16 atan(1/5)-4 atan(1/239).
# tan(4 atan(1/5)-atan(1/239))=1, with angle in (0,pi/2).
def tan_double(a):
    return 2*a/(1-a*a)
a4=tan_double(tan_double(Q(1,5)))
assert (a4-Q(1,239))/(1+a4*Q(1,239)) == 1
lo5,hi5=arctan_bounds(Q(1,5));lo239,hi239=arctan_bounds(Q(1,239))
pi_lo,pi_hi=16*lo5-4*hi239,16*hi5-4*lo239
assert Q(3)<pi_lo<pi_hi<Q(22,7)
assert Q(29,50) < Q(12,5)**2/pi_hi**2
assert Q(241,100)**2/pi_lo**2 < Q(3,5)
print('THRESHOLD: certified 29/50 < j^2/pi^2 < 3/5, consistent with 17/100 sufficient and 3/4 obstruction')

# An explicit finite-mode witness at c=2/3, beyond the proposed threshold.
# P(t)=cos(sqrt(1/3)t)+(9/25)cos(sqrt(10/3)t).
# Rational Taylor bounds on a rational grid, plus a global Lipschitz bound.
lam1,lam2,ratio = Q(1,3),Q(10,3),Q(9,25)
steps=ceil(Q(22,7)*256)
lower=min(cosine_lower(lam1,Q(i,256))+
          ratio*cosine_lower(lam2,Q(i,256)) for i in range(steps+1))
grid_lower=Q((lower*10**6).numerator//(lower*10**6).denominator,10**6)
# sqrt(1/3)<1, sqrt(10/3)<2, and |sin|<=1.
lipschitz=1+2*ratio
uniform_lower=grid_lower-lipschitz/256
assert uniform_lower > Q(1,20)
print('FINITE WITNESS: c=2/3, coefficients (1,9/25),',steps+1,'rational grid points')
print('FINITE WITNESS: grid lower >=',grid_lower,'; uniform lower >=',uniform_lower,'>1/20 on [-pi,pi]')
x0=sp.pi/sp.sqrt(2)
assert (1/sp.sqrt(2)).is_irrational and (2/sp.sqrt(2)).is_irrational
v = sp.cos(sp.sqrt(sp.Rational(1,3))*t)*sp.sin(x)/sp.sin(x0) + \
    sp.Rational(9,25)*sp.cos(sp.sqrt(sp.Rational(10,3))*t)*sp.sin(2*x)/sp.sin(2*x0)
assert sp.simplify(sp.diff(v,t,2)-sp.diff(v,x,2)-sp.Rational(2,3)*v) == 0
assert sp.simplify(v.subs(x,0)) == sp.simplify(v.subs(x,sp.pi)) == 0
assert sp.simplify(v.subs(x,x0) - (sp.cos(sp.sqrt(sp.Rational(1,3))*t)+
    sp.Rational(9,25)*sp.cos(sp.sqrt(sp.Rational(10,3))*t))) == 0
print('FINITE WITNESS: exact PDE, Dirichlet boundaries and trace at x0=pi/sqrt(2) OK')
h=sp.symbols('h',nonnegative=True)
w=sp.cosh(h*t)*sp.sin(x)
assert sp.simplify(sp.diff(w,t,2)-sp.diff(w,x,2)-(1+h*h)*w) == 0
print('HIGH-POTENTIAL FAMILY: c=1+h^2, cosh(h*t)*sin(x), exact PDE OK (h>=0)')
print('LIMIT: finite checks support the route; endpoint and all-parameter equivalence are not proved here')
print('ALL_CHECKS_OK')
