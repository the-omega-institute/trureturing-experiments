#!/usr/bin/env python3
"""Exact finite certificates and algebraic checks for the fixed-order draft.

This independently recomputes source data. No source checker is imported.
ALL_CHECKS_OK covers the stated finite checks only. No positive order-ten
witness or analytic existence theorem is certified by this program.
"""
from fractions import Fraction
from functools import reduce
from itertools import combinations, product
from math import gcd, lcm
import json
import sys

B = [
    [4,3,5,5,4,3,5,4,2,3], [4,5,5,2,2,3,4,5,3,5],
    [5,5,3,4,5,4,5,3,5,3], [2,3,4,5,2,4,4,3,3,3],
    [4,3,3,2,3,4,2,2,5,2], [5,4,5,5,2,5,5,4,3,5],
    [2,4,2,3,5,3,5,4,5,2], [4,5,2,4,3,4,2,2,5,4],
    [4,5,2,2,3,4,3,5,4,4],
]
SOURCE_C = [4725,18,-1432,4682,-7509,-7230,504,4797,5531,-5416]
SOURCE_D = [340739472931722111,28733076063275820,-369435421691654920,
            124984744924291244,-1235085791262920481,-722001405221789442,
            791954706125920512,183749257820706315,484371424027128179,
            -300611788966121134]
LOCAL = [
    (2,2,[0,0,1,0,1,1,0,1,1,0],(0,4),1),
    (3,27,[8,16,18,24,1,6,9,23,19,16],(2,8),1),
    (5,5,[3,0,2,2,2,3,3,1,2,2],(0,1),3),
    (7,7,[0,4,0,1,2,1,4,2,1,2],(1,2),5),
    (11,11,[7,2,10,5,2,7,1,7,4,8],(0,2),2),
]


def need(p, message):
    if not p:
        raise ArithmeticError(message)


def det(a):
    """Fraction-free elimination; every integer division is checked."""
    n = len(a)
    need(all(len(r) == n for r in a), 'square determinant input')
    if not n:
        return 1
    a = [list(r) for r in a]
    previous, sign = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k,n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        current = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value, remainder = divmod(current*a[i][j]-a[i][k]*a[k][j], previous)
                need(remainder == 0, 'inexact Bareiss division')
                a[i][j] = value
            a[i][k] = 0
        previous = current
    return sign*a[-1][-1]


def det_fraction(a):
    """Independent rational Gaussian elimination."""
    a = [[Fraction(x) for x in row] for row in a]
    value = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j,len(a)) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        diagonal = a[j][j]
        value *= diagonal
        for i in range(j+1,len(a)):
            multiplier = a[i][j]/diagonal
            for k in range(j+1,len(a)):
                a[i][k] -= multiplier*a[j][k]
    return value


def cofactors(a):
    n = len(a)
    return [[(-1)**(i+j)*det([r[:j]+r[j+1:] for k,r in enumerate(a) if k != i])
             for j in range(n)] for i in range(n)]


def cube(a):
    return [[x**3 for x in row] for row in a]


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def source_checks():
    c, d = [], []
    for frame,out in [(B,c),(cube(B),d)]:
        for j in range(10):
            minor = [r[:j]+r[j+1:] for r in frame]
            value = det(minor)
            need(value == det_fraction(minor), 'independent determinant algorithms')
            out.append((-1)**j*value)
    need(c == SOURCE_C and d == SOURCE_D, 'source cofactor table')
    need(reduce(gcd,c) == reduce(gcd,d) == 1, 'primitive source cofactors')
    need(-159*c[0]+41658*c[1]-c[2] == 1, 'source Bezout identity')
    gs = [gcd(d[i],d[j]) for i,j in combinations(range(10),2)]
    need(all(840 % g == 0 for g in gs) and reduce(lcm,gs) == 840, 'pair gcds')
    need(all(dot(c,r) == dot(d,[x**3 for x in r]) == 0 for r in B), 'source rows on cone')
    reduced = lambda x,i,j: c[i]*d[j]*x[j]**2-c[j]*d[i]*x[i]**2
    need(reduced(B[0],0,1) == 1123741091386468277532, 'source real Jacobian minor')
    print('FIDELITY: source Sections 2-3: all 20 cofactors, two determinant algorithms,')
    print('          primitivity, Bezout, 45 gcds, and nine homogeneous rows verified.')
    print('REAL_REDUCED_JACOBIAN_MINOR',reduced(B[0],0,1))
    for p,m,x,(i,j),expected in LOCAL:
        need(dot(c,x)%m == dot(d,[a**3 for a in x])%m == 1, 'source local congruences')
        need(reduced(x,i,j)%p == expected, 'source local reduced minor')
        if p != 3:
            need(3*expected % p != 0, 'ordinary Hensel minor')
        else:
            xi = Fraction(1-sum(c[k]*x[k] for k in range(10) if k != i),c[i])
            residual = sum(d[k]*(xi if k == i else Fraction(x[k]))**3 for k in range(10))-1
            derivative = 3*d[j]*x[j]**2-3*Fraction(d[i]*c[j],c[i])*xi**2
            need(residual.denominator%3 != 0 and residual.numerator%27 == 0, '3-adic residual')
            need(derivative.denominator%3 != 0 and derivative.numerator%3 == 0
                 and derivative.numerator%9 != 0, '3-adic derivative valuation')
        print(f'LOCAL: p={p}, modulus={m}, reduced_minor={expected}')
    return c,d


def row_relation_check(e,u,slot):
    n = len(e)+1
    y = [sum(u[k]*e[k][j] for k in range(n-1)) for j in range(n)]
    a = [list(r) for r in e]
    a.insert(slot,y)
    c = cofactors(a)
    need(det(a) == 0, 'inserted row lies in row span')
    for k in range(n-1):
        index = k if k < slot else k+1
        need(c[index] == [-u[k]*z for z in c[slot]], 'signed cofactor proportionality')


def identity_checks():
    count = 0
    for n in [2,3]:
        for entries in product([1,2],repeat=n*(n-1)):
            e = [list(entries[k*n:(k+1)*n]) for k in range(n-1)]
            for u in product([-1,1,2],repeat=n-1):
                for slot in range(n):
                    row_relation_check(e,u,slot)
                    row_relation_check(cube(e),u,slot)
                    count += 2
    for slot in range(10):
        row_relation_check(B,list(range(1,10)),slot)
        row_relation_check(cube(B),list(range(1,10)),slot)
        count += 2
    print('COFACTOR_IDENTITY_CASES',count)
    print('SCOPE: every 1x2 and 2x3 frame with entries in {1,2}; every relation')
    print('       coefficient in {-1,1,2}; every insertion slot; and the source 9x10 frame.')
    base = [[1]+[0]*9]+B
    originals = [cofactors(base),cofactors(cube(base))]
    scaling_count = 0
    for slot in range(10):
        for scale in [2,10,100]:
            a = [r[:] for r in base]
            a[slot] = [scale*x for x in a[slot]]
            for power,matrix in [(1,a),(3,cube(a))]:
                cs = cofactors(matrix)
                for k in range(10):
                    expected = 1 if k == slot else scale**power
                    need(cs[k] == [expected*v for v in originals[power == 3][k]],
                         'row replacement exact cofactor homogeneity')
                scaling_count += 1
    print('ORDER_TEN_SCALING_CASES',scaling_count,'scales=2,10,100')
    from sympy import symbols, expand
    a = symbols('a0:3'); c = symbols('c0:3'); d = symbols('d0:3'); lam = symbols('lam')
    expression = sum(a[j]*(3*d[j]*a[j]**2-lam*c[j]) for j in range(3))
    need(expand(expression-(3*sum(d[j]*a[j]**3 for j in range(3))
                            -lam*sum(c[j]*a[j] for j in range(3)))) == 0,
         'Jacobian singularity scalar identity')
    print('SINGULARITY_IDENTITY: sum a_j(3 d_j a_j^2 - lambda c_j) = 3 C - lambda L.')


def witnesses(path=None):
    if path is None:
        for bound in [2,10,100]:
            print(f'REQUIRED_WITNESS_NOT_SUPPLIED n=10 K={bound}')
        print('No determinant-one positive order-ten witness is claimed.')
        return
    rows = json.load(open(path))
    for bound in [2,10,100]:
        a = rows[str(bound)]
        need(len(a) == 10 and all(len(r) == 10 for r in a), 'witness order')
        need(all(type(x) is int and x > max(bound,1) for r in a for x in r), 'witness bound')
        ds = det(a),det(cube(a))
        need(ds == (1,1), 'witness determinants')
        need(det_fraction(a) == det_fraction(cube(a)) == 1, 'independent witness determinants')
        print('WITNESS n=10 K=',bound,'minimum=',min(map(min,a)),'determinants=',ds)


def main():
    source_checks()
    identity_checks()
    witnesses(sys.argv[1] if len(sys.argv) > 1 else None)
    print('LIMITATION: finite identities do not establish the quantified existence claim.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
