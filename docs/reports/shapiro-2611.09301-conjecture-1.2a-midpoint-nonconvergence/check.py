#!/usr/bin/env python3
"""Exact checks for Shapiro arXiv:2611.09301v1, Conjecture 1.2(a).

All mathematical comparisons use rational/integer arithmetic or symbolic algebra.
Finite checks support a proposed general argument; they do not certify its Lean proof.
"""
from fractions import Fraction as F
from math import factorial, comb
from functools import reduce
from operator import mul
import sys
import sympy as S

assert sys.version_info >= (3, 12)
t = S.Symbol('t')
EPS = F(1, 559872)


def prod(xs):
    return reduce(mul, xs, F(1))


def coefficient(k, x):
    if k == 0:
        return 1 - 2*x
    if k == 1:
        return x - F(1, 4)
    if k == -1:
        return x - F(3, 4)
    return F(0)


def matrix(n):
    return S.Matrix(n, n, lambda i, j:
                    S.Rational(coefficient(j-i, F(i+j, 2*n))))


def charpoly(n):
    """Continuant built directly from the literal midpoint entries."""
    prev, cur = S.Poly(1, t, domain=S.QQ), S.Poly(t-1, t, domain=S.QQ)
    if n == 0:
        return prev
    for i in range(1, n):
        diagonal = coefficient(0, F(i, n))
        upper = coefficient(1, F(2*i-1, 2*n))
        lower = coefficient(-1, F(2*i-1, 2*n))
        prev, cur = cur, S.Poly(t-diagonal, t)*cur - upper*lower*prev
    return cur


def rising(w, n):
    return prod(w+j for j in range(n))


def q(r):
    return prod(F(2*j-1, 2) for j in range(1, r+1))


def zero_order(p):
    return min(monomial[0] for monomial, c in p.terms() if c != 0)


print('Exact arithmetic: fractions.Fraction and sympy', S.__version__)
print('Source: Shapiro, arXiv:2611.09301v1, Conjecture 1.2(a), equation (1.1)')
print('Sampling: a[j-i]((i+j)/(2*N)), 0 <= i,j < N')
# Fidelity: source Section 1.3 gives b+2*c*cos(j*pi/(m+1)).
b, c = S.symbols('b c')
frozen = S.Matrix([[b,c,0],[c,b,c],[0,c,b]])
source_roots = [b+2*c*S.cos(j*S.pi/4) for j in range(1,4)]
assert S.expand(frozen.charpoly(t).as_expr() - prod(t-z for z in source_roots)) == 0
print('FIDELITY_OK: symmetric frozen 3-by-3 charpoly = (t-b)*((t-b)^2-2*c^2)')

polys = {n: charpoly(n) for n in range(1, 51)}
for n in range(1, 13):
    assert S.expand(matrix(n).charpoly(t).as_expr() - polys[n].as_expr()) == 0
print('DENSE_CROSSCHECK_OK: literal matrices vs independent continuants, N=1..12')

# Check the shifted recurrence formula for every leading size, including zeros.
checks = 0
for n in range(1, 51):
    u, v = F(1,2)-F(n,4), F(1,2)-F(3*n,4)
    pprev, pcur = F(1), u+v-1
    sequence = [pprev, pcur]
    for j in range(1,n):
        pprev, pcur = pcur, (2*j+u+v-1)*pcur-(j+u-1)*(j+v-1)*pprev
        sequence.append(pcur)
    for j in range(n+1):
        formula = ((u-1)*rising(u,j)-(v-1)*rising(v,j))/(u-v)
        assert sequence[j] == formula
        checks += 1
    assert sequence[n] == F(polys[n].nth(0))*n**n
print('SHIFTED_CONTINUANT_OK:', checks, 'leading determinants, N=1..50')

for r in range(1, 129):
    assert q(r) >= F(factorial(r),2**r)
    assert factorial(r) >= F(r,3)**r
    assert sum((F(comb(r,j),r**j) for j in range(r+1)),F(0)) <= 3
print('FACTORIAL_BOUNDS_OK: r=1..128; epsilon =', EPS)

for n in range(1,51):
    for i in range(n):
        row = abs(coefficient(0,F(i,n)))
        if i+1<n:
            row += abs(coefficient(1,F(2*i+1,2*n)))
        if i>0:
            row += abs(coefficient(-1,F(2*i-1,2*n)))
        assert row <= F(5,2)
    for i in range(n-1):
        x=F(2*i+1,2*n)
        assert coefficient(1,x)-coefficient(-1,x) == F(1,2)
        assert coefficient(1,x) != 0 or coefficient(-1,x) != 0
print('ROW_BOUND_AND_CONNECTED_PATH_OK: N=1..50')

print('SUBSEQUENCE_CERTIFICATES: f(z)=max(0,1-|z|/epsilon)')
for m in range(1,13):
    n = 4*m
    p = polys[n]
    determinant = F(p.nth(0))  # N is even.
    expected = (-1)**m*q(m)*q(3*m)/n**n
    assert determinant == expected
    assert abs(determinant) >= F(1,12**n)
    assert all(coefficient(k,F(2*i+1,2*n)) != 0
               for i in range(n-1) for k in [-1,1])
    # If more than N/4 roots were in the closed epsilon disk, their product
    # and the radius-3 row bound would contradict this exact strict inequality.
    forbidden = n//4+1
    assert EPS**forbidden*3**(n-forbidden) < F(1,12**n)
    # A stronger finite check: the constant term dominates on the closed disk.
    tail = sum((abs(F(p.nth(j)))*EPS**j for j in range(1,n+1)),F(0))
    no_small_roots = tail < abs(F(p.nth(0)))
    assert no_small_roots

    nn = 4*m+2
    qq = 2*m+1
    pn = polys[nn]
    mat = matrix(nn)
    assert mat[m,m+1] == 0
    assert mat[3*m+2,3*m+1] == 0
    K = (nn*mat)[m+1:3*m+2,m+1:3*m+2]
    direct = S.zeros(qq)
    for j in range(qq):
        direct[j,j] = qq-1-2*j
        if j+1<qq:
            direct[j,j+1] = j+1
            direct[j+1,j] = j-(qq-1)
    assert K == direct
    translation = S.Matrix(qq,qq,lambda i,j:
                           comb(j,i)*(-1)**(j-i) if i<=j else 0)
    shift = S.zeros(qq)
    for j in range(qq-1):
        shift[j+1,j] = j-qq+1
    assert translation.det() == 1
    assert K*translation == translation*shift
    assert K.charpoly(t).as_expr() == t**qq
    assert zero_order(pn) >= qq
    assert pn.rem(S.Poly(t**qq,t)).is_zero
    print('m=%d N=%d det=%s roots_in_epsilon_disk=0 integral_f=0; '
          'N=%d zero_multiplicity=%d integral_f>=1/2' %
          (m,n,determinant,nn,zero_order(pn)))

for n in [4,6,8,10,12,14]:
    print('EXACT_CHARPOLY N=%d:' % n, S.factor(polys[n].as_expr()))
assert 8*matrix(4) == S.Matrix([[8,-1,0,0],[-5,4,1,0],[0,-3,0,3],[0,0,-1,-4]])
assert polys[4].nth(0) == -S.Rational(15,4096)
assert S.expand(polys[6].as_expr() - t**3*(t+S.Rational(2,3))*
                (t*t-S.Rational(5,3)*t+S.Rational(5,9))) == 0
print('SMALL_INSTANCES_OK: matrices N=4 and exact spectrum N=6')
for n in range(1,9):
    for rr, ss in [(S.Rational(2),S.Rational(3)),(S.I,1+S.I),(1+S.I,S.Rational(-2,3))]:
        transformed = ss*S.eye(n)+rr*matrix(n)
        expected = S.expand(rr**n*polys[n].as_expr().subs(t,(t-ss)/rr))
        assert S.expand(transformed.charpoly(t).as_expr()-expected) == 0
print('AFFINE_FAMILY_OK: three exact nonzero complex scales, N=1..8')
print('Scope: finite checks; the proposed infinite-family argument is not a Lean proof.')
print('ALL_CHECKS_OK')
