#!/usr/bin/env python3
"""Exact supporting checks for the three-state Bargmann-density proposal.

No numerical quadrature or floating point is used. Finite checks do not prove
Haar disintegration, dominated convergence, or the two proposed asymptotics.
"""
from fractions import Fraction as F
from math import comb, factorial, prod
import sympy as S


def rising(a, k):
    return prod(range(a, a + k))


def beta_integer(a, b):
    return F(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))


def source_moment(d, m, n):
    """Zhang Theorem 6.7, TeX label eq:mixed, with s=m+n."""
    total = sum(comb(m, j) * comb(n, j) * factorial(m + n - j)**2
                * factorial(j) * rising(d - 1, j)
                for j in range(min(m, n) + 1))
    return F(total, rising(d, m + n)**2)


def disk_moment(d, m, n):
    """Integrate the proposed conditional disk law Z=uv+R*xi exactly."""
    answer = F(0)
    for j in range(min(m, n) + 1):
        # Independent u,v have Beta(1,d-1) law. xi is one coordinate
        # of an independent Haar vector in C^(d-1).
        uv = (d - 1) * beta_integer(m + n - j + 1, d + j - 1)
        radial = F(factorial(j) * factorial(d - 2), factorial(d + j - 2))
        answer += comb(m, j) * comb(n, j) * uv**2 * radial
    return answer


u, v, x, y = S.symbols('u v x y', positive=True)
r2 = x*x + y*y
R2 = u*(1-u)*v*(1-v)
D = 1 + 2*x - u - v - r2/(u*v)
assert S.factor((1-u)*(1-v) - ((x-u*v)**2+y*y)/(u*v) - D) == 0
print('DENSITY: conditional disk deficit identity verified symbolically')

# Source fidelity: evaluate the source mixed-moment formula in the qubit case.
ez = source_moment(2, 1, 0)
ez2 = source_moment(2, 2, 0)
absz2 = source_moment(2, 1, 1)
assert (ez, (absz2+ez2)/2, (absz2-ez2)/2) == (F(1,4), F(1,8), F(1,72))
print('FIDELITY: d=2 E[X]=1/4, E[X^2]=1/8, E[Y^2]=1/72 (source elementary marginals)')

cases = 0
for d in range(3, 10):
    p = d-3
    # Integrate in z first: the conditional disk contributes pi*uv
    # *((1-u)*(1-v))^(p+1)/(p+1). pi cancels C_d.
    mass = F((d-1)**2*(d-2), (p+1)*(p+2)**2)
    assert mass == 1
    for m in range(5):
        for n in range(5):
            assert disk_moment(d,m,n) == source_moment(d,m,n)
            cases += 1
print(f'DENSITY: normalization for d=3..9 and all {cases} mixed moments with m,n=0..4 agree exactly')

# The logarithmic coordinate map has measure ds dt and a triangular domain.
L, h = S.symbols('L h', positive=True)
c = S.log(1/h)
for p in range(9):
    upper = (1+2*S.exp(-L))**p * (2*L+S.log(1+2*S.exp(-L)))**2/(2*L**2)
    lower = (1-2*S.exp(-L)-3*h)**p * (2*L-3*c)**2/(2*L**2)
    assert S.limit(upper,L,S.oo) == 2
    lower_limit = S.limit(lower,L,S.oo)
    assert S.simplify(lower_limit-2*(1-3*h)**p) == 0
    assert S.limit(lower_limit,h,0,dir='+') == 2
print('ORIGIN: exact upper/lower triangle limits both give 2 for p=0..8; angular bounds use |Re z|<=|z|')

s, t, q = S.symbols('s t q', real=True)
Phi = S.exp(s)+S.exp(t)+S.exp(-s-t)-3
Q = s*s+s*t+t*t
assert S.expand(S.series(Phi.subs({s:q*s,t:q*t}),q,0,3).removeO()-q*q*Q) == 0
assert S.hessian(Phi,(s,t)).subs({s:0,t:0}) == S.Matrix([[2,1],[1,2]])
alpha,beta = s+t/2,S.sqrt(3)*t/2
assert S.expand(alpha**2+beta**2-Q) == 0
assert S.Matrix([alpha,beta]).jacobian([s,t]).det() == S.sqrt(3)/2
rho = S.symbols('rho', nonnegative=True)
for p in range(9):
    radial_integral = S.integrate(rho*(1-rho*rho)**p,(rho,0,1))
    assert radial_integral == S.Rational(1,2*(p+1))
    d=p+3
    # C_d * a^p * [2*pi/(sqrt(3)*(p+1))] * (delta/a)^(p+1).
    coefficient = S.cancel(S.Rational((d-1)**2*(d-2),1)/S.pi * 2*S.pi/(S.sqrt(3)*(p+1)))
    assert S.simplify(coefficient-2*(d-1)**2/S.sqrt(3)) == 0
print('BOUNDARY: Hessian, quadratic Jacobian, beta integrals, and proposed coefficient verified exactly for d=3..11')

# Rational boundary points: r=q^3, a=q^2, x=(3q^2-1)/2.
count=0
for den in range(2,17):
    for num in range((den+1)//2,den+1):
        a=F(num,den)**2
        r=F(num,den)**3
        xb=(3*a-1)/2
        yb2=r*r-xb*xb
        assert yb2>=0
        assert 27*(xb*xb+yb2)==(1+2*xb)**3
        assert -F(1,8)<=xb<=1
        assert r>=F(1,8)
        assert a>=F(1,4)
        count+=1
assert (F(1,2)**3, (3*F(1,2)**2-1)/2)==(F(1,8),-F(1,8))
assert (F(1)**3,(3*F(1)**2-1)/2)==(F(1),F(1))
eps=S.symbols('eps',positive=True)
assert S.limit((1+2*(1-eps)-3*(1-eps)**S.Rational(2,3))/eps**2,eps,0,dir='+')==S.Rational(1,3)
print(f'ENDPOINTS: {count} exact rational boundary parameter checks, including -1/8 and 1; cusp delta/(1-x)^2 -> 1/3')
print('SCOPE: exact supporting identities and limits; no Lean proof and no full certification of the density asymptotics')
print('ALL_CHECKS_OK')
