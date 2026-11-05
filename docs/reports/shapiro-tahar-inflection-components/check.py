#!/usr/bin/env python3
"""Exact algebra and uniform boundary certificates; no sampled plots or floats.

The topological implications of these certificates are explained in README.md.
This program does not decide connected components of arbitrary real curves.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import sympy as S

z = S.Symbol('z')
x, y = S.symbols('x y', real=True)
I = S.I


def exact_zero(e):
    assert S.cancel(S.expand(e)) == 0, e


def cleared_imaginary(r):
    q, p = map(S.expand, S.fraction(S.cancel(r)))
    w = S.diff(q, z) * p - q * S.diff(p, z)
    pc = S.expand(p.subs(z, x + I*y))
    wc = S.expand(w.subs(z, x + I*y))
    return S.expand((wc*S.conjugate(pc)**2-S.conjugate(wc)*pc**2)/(2*I))


# Source Section 3: W/P^2 is the derivative of Q/P.
for q in [1+z, z**2+I*z+1, z**4+2*I]:
    for p in [z, z**2+1, z**3-z+I]:
        exact_zero(S.diff(q/p, z)-(S.diff(q,z)*p-q*S.diff(p,z))/p**2)
print('FIDELITY: Section 3 quotient-derivative identity: 9 exact polynomial pairs')

# Source Section 11: R=1/z and R=z+1/z give the two coordinate axes.
# 2*x*y=0 is a connected unbounded union, including the pole at 0.
for r in [1/z, z+1/z]:
    exact_zero(cleared_imaginary(r)-2*x*y)
# Source Section 11: R=z^d gives Im(z^(d-1))=0, d-1 lines through 0.
for d in range(2,9):
    exact_zero(cleared_imaginary(z**d)-d*S.im((x+I*y)**(d-1)).expand(complex=True))
# Source Theorem 4.1: derivative has exact leading pole order s+1.
for order in range(1,6):
    a = 2+3*I
    exact_zero(S.cancel(z**(order+1)*S.diff(a/z**order+z+I,z))+order*a-(1)*z**(order+1))
print('FIDELITY: Section 11 pole/noncompact examples; polynomial degrees 2..8; local orders 1..5')


radius = F(1,32)

def perturbation_bound(poles, orders, inserted, j):
    return sum((F(orders[k],1)/(radius if k==j else F(abs(poles[k]-poles[j]))-radius)**(orders[k]+1)
                for k in inserted), F(0))


def polynomial_margin(p, degree):
    a = F(p+1)
    k = degree-1
    deviation = sum((F(comb(k,l))*a**(k-l)*radius**l for l in range(1,k+1)),F(0))
    return degree*(a**k-deviation)


def proper_margin(p, order):
    a=F(p);k=order+1
    deviation=((a+radius)**k-a**k)/(a**k*(a-radius)**k)
    return order*(1/a**k-deviation)


def certify(poles, orders, infinity_order):
    n=len(poles)
    if infinity_order:
        r0=I*(z+1)**infinity_order
        inserted=list(range(n))
        margins=[polynomial_margin(poles[j],infinity_order) for j in inserted]
    else:
        assert n>=1
        r0=I/z**orders[0]
        inserted=list(range(1,n))
        margins=[proper_margin(poles[j],orders[0]) for j in inserted]
    bounds=[perturbation_bound(poles,orders,inserted,j) for j in inserted]
    assert all(a>0 for a in margins)
    eps=min((a/(2*(b+1)) for a,b in zip(margins,bounds)),default=F(1))
    assert eps>0
    assert all(a-eps*b>0 for a,b in zip(margins,bounds))
    r=r0+S.Rational(eps.numerator,eps.denominator)*sum((1/(z-poles[k])**orders[k] for k in inserted),S.Integer(0))
    q,p=map(lambda v:S.Poly(v,z,extension=I),S.fraction(S.cancel(r)))
    assert not q.as_expr().has(S.Float) and not p.as_expr().has(S.Float)
    assert max(q.degree()-p.degree(),0)==infinity_order
    expected_p=S.prod((z-a)**s for a,s in zip(poles,orders))
    exact_zero(p.monic().as_expr()-expected_p)
    assert all(q.eval(a)!=0 for a in poles)
    # This equals the desired maximum once the prose topological upper bound
    # and boundary isolation lemma are proved; it is a certificate count here.
    assert len(inserted)==(n if infinity_order else n-1)
    return len(inserted)

cases=certificates=0
for n in range(4):
    for orders in product(range(1,4),repeat=n):
        for mi in range(4):
            if n==0 and mi==0:continue
            certificates+=certify(tuple(range(n)),orders,mi)
            cases+=1
print(f'SETTLEMENT CHECK: {cases} prescribed-divisor constructions, {certificates} uniform disk certificates')
print('Domain checked: n=0..3, each finite pole order=1..3, infinity order=0..3, nonzero divisor; positions=0,..,n-1')

# Candidate's first explicit check, with an intentionally coarse bound.
pert=F(16,100)+2*F(16,9*100)
assert pert==F(44,225) and 1-pert==F(181,225)>0
r=I*z+S.Rational(1,100)*(1/z+1/(z-1)+1/(z-2))
q,p=S.fraction(S.cancel(r))
assert S.degree(q,z)-S.degree(p,z)==1
assert all(S.expand(q).subs(z,a)!=0 for a in [0,1,2])
print('EXPLICIT polynomial background: perturbation <=44/225; Im derivative >=181/225 on three radius-1/4 boundaries')

# On |z-1|=1/4: x>=3/4, |y|<=1/4 and |z|<=5/4.
# Re(1/z^2)=(x^2-y^2)/(x^2+y^2)^2 >=(1/2)/(5/4)^4.
lower=(F(3,4)**2-F(1,4)**2)/F(5,4)**4
pert=F(1,100)/F(1,4)**2
assert lower==F(128,625) and pert==F(100,625)
assert lower-pert==F(28,625)>0
r=I/z+S.Rational(1,100)/(z-1)
q,p=S.fraction(S.cancel(r))
assert S.degree(q,z)<=S.degree(p,z)
assert S.expand(q).subs(z,0)!=0 and S.expand(q).subs(z,1)!=0
print('EXPLICIT proper background: Im derivative <=-28/625 on radius-1/4 boundary about 1')
print('LIMIT: exact algebra/boundary checks support the construction; general topology remains a future Lean obligation')
print('ALL_CHECKS_OK')
