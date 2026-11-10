#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact source fidelity and K_{m,m} curvature/diversity-limit checks.

Requires Python >= 3.12 and SymPy. All evidence uses integer/rational
arithmetic or symbolic polynomial identities; no floating point is used.
"""
from fractions import Fraction
from itertools import product
import sys
import sympy as S

assert sys.version_info >= (3, 12)


def zero(expression):
    assert S.cancel(S.expand(expression)) == 0, expression


def distance(m):
    return S.Matrix(2*m, 2*m, lambda i, j:
                    0 if i == j else 2 if i//m == j//m else 1)


def source_fidelity():
    # Roff--Willerton v2, the five-point example (TeX label eg:five-points).
    a, b = S.symbols('a b')
    D = S.Matrix([[0,a,a,1,1], [a,0,a,1,1], [a,a,0,1,1],
                  [1,1,1,0,b], [1,1,1,b,0]])
    zero(D.det() - 2*(3-a*b)*a*a*b)
    den = 4*a + 3*b - 12
    v = S.Matrix([b-2,b-2,b-2,2*a-3,2*a-3])/den
    c = 2*(a*b-3)/den
    zero(sum(v)-1)
    for entry in D*v-S.ones(5,1)*c:
        zero(entry)
    exceptional = D.subs({a:S.Rational(3,2), b:2})
    assert exceptional*(S.ones(5,1)/5) == S.ones(5,1)
    print('FIDELITY: source five-point determinant, generic gauging, exceptional gauging OK')


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first, length-1):
                yield (first,)+rest


def family_case(m):
    D = distance(m)
    n = 2*m
    # The actual finite metric, checked independently of spectral algebra.
    for i,j,k in product(range(n), repeat=3):
        assert D[i,j] >= 0 and D[i,j] == D[j,i]
        assert (D[i,j] == 0) == (i == j)
        assert D[i,k] <= D[i,j]+D[j,k]
    lam = S.symbols('lambda')
    expected = (lam-(3*m-2))*(lam-(m-2))*(lam+2)**(2*m-2)
    zero(D.charpoly(lam).as_expr()-expected)
    assert D.det() != 0
    one = S.ones(n,1)
    u = one/n
    assert D*u == one*S.Rational(3*m-2,2*m)
    # Solve the full mass-one/constant-potential system, not just D u = D u.
    v = S.symbols('v:'+str(n))
    c = S.symbols('c')
    solution = S.linsolve(list(D*S.Matrix(v)-c*one)+[sum(v)-1], (*v,c))
    assert solution == S.FiniteSet((*([S.Rational(1,n)]*n),S.Rational(3*m-2,n)))
    # The mass-zero split direction proves this family is outside negative type.
    h = S.Matrix([1]*m+[-1]*m)
    assert sum(h) == 0
    assert (h.T*D*h)[0] == 2*m*(m-2) > 0

    r = S.symbols('r')
    p = S.symbols('p:'+str(n))
    a, b = sum(p[:m]), sum(p[m:])
    Z = D.applyfunc(lambda d:r**d)
    E = (S.Matrix(p).T*Z*S.Matrix(p))[0]
    raw = (1-r*r)*sum(x*x for x in p)+r*r*(a*a+b*b)+2*r*a*b
    zero(E-raw)
    cr = (1+(m-1)*r*r)/m
    variance = sum((x-a/m)**2 for x in p[:m])+sum((x-b/m)**2 for x in p[m:])
    decomposed = cr+2*(r-cr)*a*(1-a)+(1-r*r)*variance
    # Enforce exactly mass one, leaving every other coordinate symbolic.
    zero((E-decomposed).subs(p[-1],1-sum(p[:-1])))
    zero(r-cr-(1-r)*((m-1)*r-1)/m)
    # All signs are certified by coefficients after mapping each open interval
    # rationally from s>0. Positive coefficients imply strict positivity there.
    s = S.symbols('s',positive=True)
    low = S.Rational(1,m-1)
    r_small_scale = low+(1-low)*s/(1+s)
    r_large_scale = low*s/(1+s)
    for expr in [(r-cr).subs(r,r_small_scale),
                 (cr-r).subs(r,r_large_scale),
                 (1-r*r).subs(r,r_small_scale),
                 (1-r*r).subs(r,r_large_scale)]:
        numerator, denominator = S.fraction(S.cancel(expr))
        for poly in [S.Poly(numerator,s),S.Poly(denominator,s)]:
            assert all(coef >= 0 for coef in poly.all_coeffs())
            assert any(coef > 0 for coef in poly.all_coeffs())
    assert S.cancel((r-cr).subs(r,low)) == 0
    pA = S.Matrix([S.Rational(1,m)]*m+[0]*m)
    pB = S.Matrix([0]*m+[S.Rational(1,m)]*m)
    zero((pA.T*Z*pA)[0]-cr)
    zero((pB.T*Z*pB)[0]-cr)
    assert sum(abs(x-y) for x,y in zip(pA,u))/2 == S.Rational(1,2)
    assert sum(abs(x-y) for x,y in zip(pB,u))/2 == S.Rational(1,2)
    mix = S.Matrix([a/m]*m+[(1-a)/m]*m)
    # At the threshold all separately uniform mixtures have exactly c_r energy.
    zero(((mix.T*Z*mix)[0]-cr).subs(r,low))
    # Independently evaluate the objective on every denominator-four simplex
    # point at rational kernels on either side and exactly at the threshold.
    count = 0
    for ints in compositions(4,n):
        pp = [Fraction(x,4) for x in ints]
        aa = sum(pp[:m]); bb = 1-aa
        vv = sum((x-aa/m)**2 for x in pp[:m])+sum((x-bb/m)**2 for x in pp[m:])
        for rr in [(Fraction(1,1)+Fraction(1,m-1))/2,
                   Fraction(1,m-1),Fraction(1,2*(m-1))]:
            energy = sum(pp[i]*pp[j]*rr**int(D[i,j]) for i in range(n) for j in range(n))
            cc = (1+(m-1)*rr*rr)/m
            assert energy == cc+2*(rr-cc)*aa*bb+(1-rr*rr)*vv
            optimum = cc if rr >= Fraction(1,m-1) else (cc+rr)/2
            assert energy >= optimum
        count += 1
    print(f'FAMILY m={m}: metric, spectrum, unique positive curvature, failure of negative type,')
    print(f'  exact energy/variance identities, three scale regimes, TV=1/2; {count} rational simplex points OK')


if __name__ == '__main__':
    source_fidelity()
    for m in [3,4,5]:
        family_case(m)
    # Symbolic factorization with m arbitrary complements the fixed-size checks.
    m,r = S.symbols('m r')
    zero(r-(1+(m-1)*r*r)/m-(1-r)*((m-1)*r-1)/m)
    print('GENERAL: arbitrary-m transition factorization OK')
    print('ALL_CHECKS_OK')
