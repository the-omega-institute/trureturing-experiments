"""有限三体基态与平方占据插值在同一 Regge 形变上的一、二阶响应。"""
from fractions import Fraction
from itertools import combinations
from math import comb, isqrt

import sympy as sp
from flint import arb, ctx


def angle_jets():
    t = sp.Symbol("t")
    def squared_edge(i,j):
        width = abs(j-i)
        return sp.Rational(width*(5-width),5)+int(width in (1,4))*t
    gram = sp.Matrix(4,4,lambda i,j: (squared_edge(0,i+1)+squared_edge(0,j+1)
                                       -squared_edge(i+1,j+1))/2)
    normals = (-sp.ones(1,4)).col_join(sp.eye(4))
    products = normals*gram.inv()*normals.T
    jets = {4:[0,0],6:[0,0]}
    for i,j in combinations(range(5),2):
        cosine = -products[i,j]/sp.sqrt(products[i,i]*products[j,j])
        at_zero = cosine.subs(t,0)
        sine = sp.sqrt(1-at_zero**2)
        first = cosine.diff(t).subs(t,0)/sine
        second = cosine.diff(t,2).subs(t,0)/sine+at_zero*first**2/sine
        pair = jets[4 if at_zero==0 else 6]
        pair[0] += 24*first
        pair[1] += 24*second
    assert jets == {4:[90,-90],6:[-60*sp.sqrt(3),90*sp.sqrt(3)]}
    for indices,expected in [((0,1,3),sp.Rational(1,5)),((0,1,2),sp.Rational(3,10))]:
        i,j,k = indices
        a,b,c = squared_edge(i,j),squared_edge(j,k),squared_edge(i,k)
        area_squared = (2*(a*b+b*c+c*a)-a*a-b*b-c*c)/16
        assert area_squared.diff(t).subs(t,0)==expected


def mean_energy(n,k,x):
    return Fraction(2,3)+Fraction(comb(k,3),3*comb(n,3))-x*Fraction(k,n)


def cold(n,x):
    target = x*(n-1)*(n-2)
    k = (1+isqrt(1+4*(target.numerator//target.denominator)))//2
    if k*(k-1)<target:
        k += 1
    assert (k-1)*(k-2)<target<k*(k-1)
    assert mean_energy(n,k,x)<min(mean_energy(n,k-1,x),mean_energy(n,k+1,x))
    return Fraction(k,n),Fraction(0)


def interpolate(n,x):
    squared = n*n*x
    k = isqrt(squared.numerator//squared.denominator)
    assert k*k<squared<(k+1)**2
    weight = (squared-k*k)/(2*k+1)
    assert (1-weight)*k*k+weight*(k+1)**2 == squared
    return (k+weight)/n,Fraction(n,2*k+1)


def ball(z):
    return arb(z.numerator)/z.denominator


def shape_response(method,n):
    p4,v4 = map(ball,method(n,Fraction(1,5)))
    p6,v6 = map(ball,method(n,Fraction(3,20)))
    root = arb(3).sqrt()
    force = 90*p4-60*root*p6
    hessian = -90*p4+90*root*p6+2*(90*v4/5-60*root*v6*3/10)
    return force,hessian


def verify():
    ctx.prec = 192
    angle_jets()
    n = 1000
    target = -45/arb(5).sqrt()
    for method in (cold,interpolate):
        force,hessian = shape_response(method,n)
        print("有限槽位、规则、形状力、二阶系数、相对 Regge 二阶系数",
              n,method.__name__,force.str(12),hessian.str(12),(hessian/target).str(12))
        if method is cold:
            assert -1.02<hessian/target<-.98
        else:
            assert abs(force)<arb(".0001") and abs(hessian/target-1)<arb(".01")
    for size in (20,1000,20000):
        for x in (Fraction(1,5),Fraction(3,20)):
            p,slope = interpolate(size,x)
            r = ball(x).sqrt()
            assert abs(ball(cold(size,x)[0])-r)<arb(2)/size
            assert 0 < r-ball(p) <= 1/(4*size**2*(r+ball(p)))
            k = isqrt((size**2*x).numerator//(size**2*x).denominator)
            assert abs(2*r*ball(slope)-1)<arb(1)/(2*k+1)


if __name__ == "__main__":
    verify()
