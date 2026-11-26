"""两档有限占据的弯曲驻点、角亏与每胞平均标量曲率。"""
from fractions import Fraction
from itertools import combinations
from math import isqrt
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import arb, ctx


def simplex():
    t = sp.Symbol("t",positive=True)
    def squared_edge(i,j):
        if i==j:
            return 0
        return t if abs(i-j) in (1,4) else sp.S.One
    gram = sp.Matrix(4,4,lambda i,j: (squared_edge(0,i+1)+squared_edge(0,j+1)
                                       -squared_edge(i+1,j+1))/2)
    normals = (-sp.ones(1,4)).col_join(sp.eye(4))
    form = normals*gram.inv()*normals.T
    assert sp.factor(gram.det()) == 5*(t*t-3*t+1)**2/16
    for i,j in combinations(range(5),2):
        assert sp.factor(form[i,i]-form[j,j]) == 0
        expected = (3-2*t if j-i in (1,4) else 3*t-2)/(2*(1+t))
        assert sp.factor(-form[i,j]/form[i,i]-expected)==0


def readings(k4,k6,n):
    ratio = arb(k4)/k6
    t = 2*(1-ratio**2)+(4*(1-ratio**2)**2+ratio**2).sqrt()
    b = 4*arb(k6)/(n*(4*t-1).sqrt())
    assert arb(".64")<t<arb(".69")
    for k in (k4,k6):
        x = Fraction(k*k,n*n)*(n-1)*(n-2)
        assert (k-1)*(k-2)<x<k*(k-1)
    four = 2*arb.pi()-4*((3*t-2)/(2*(1+t))).acos()
    six = 2*arb.pi()-6*((3-2*t)/(2*(1+t))).acos()
    volume = arb(5).sqrt()*b*b*(3*t-t*t-1)/4
    R = (30*k4*four+20*k6*six)/n
    curvature = 2*R/volume
    assert four*six<0 and curvature<0
    return t,b,four,six,curvature


def verify():
    ctx.prec = 224
    simplex()
    cold = run_path(str(Path(__file__).with_name("ch39-hinge-occupation-response.py")))["cold"]
    t = sp.Symbol("t")
    flat = sp.Rational(2,3)
    ratio_squared = t*(4-t)/(4*t-1)
    # 2t²−t+2 = 2(t−1/4)²+15/8，故面积平方比严格递减。
    assert sp.factor((4*t-1)**2 * ratio_squared.diff(t) + 2*(2*t**2-t+2)) == 0
    slope = ratio_squared.diff(t).subs(t, flat)
    assert slope == -sp.Rational(8,5)
    cos_four = (3*t-2)/(2*(1+t))
    cos_six = (3-2*t)/(2*(1+t))
    assert sp.factor(cos_six.diff(t)) == -5/(2*(1+t)**2)
    area_four = sp.sqrt(t*(4-t))/4
    area_six = sp.sqrt(4*t-1)/4
    action = (30*area_four*(2*sp.pi-4*sp.acos(cos_four))
              + 20*area_six*(2*sp.pi-6*sp.acos(cos_six)))
    assert sp.simplify(action.subs(t,flat)) == 0
    assert sp.simplify(action.diff(t).subs(t,flat)) == 0
    quadratic = sp.simplify(action.diff(t,2).subs(t,flat)/2)
    assert quadratic == -27/sp.sqrt(5)
    unit_volume = sp.sqrt(5)*(3*flat-flat**2-1)/4
    # 固定体积使 A_i² 正比于 h_i/v，以下两式给相应形状导数的负分子。
    volume_polynomial = 3*t-t**2-1
    for area_polynomial, numerator in [(t*(4-t), 4-2*t-t**2),
                                       (4*t-1, 1+2*t-4*t**2)]:
        assert sp.expand(volume_polynomial.diff(t)*area_polynomial
                         - volume_polynomial*area_polynomial.diff(t) - numerator) == 0
    flat_b = sp.Rational(6,5)
    flat_p6 = flat_b*area_six.subs(t,flat)
    limiting_coefficient = sp.simplify(-2*quadratic / (unit_volume*flat_b)
                                       / (3*slope)**2 / flat_p6**4)
    assert limiting_coefficient == 125
    for n in (1000,10000,100000):
        p4,p6 = [cold(n,x)[0] for x in (Fraction(1,5),Fraction(3,20))]
        k4,k6 = int(n*p4),int(n*p6)
        row = readings(k4,k6,n)
        if n == 1000:
            prep_force = -2*(30*p4**3+20*p6**3)
            assert prep_force == -Fraction(604013,78125)
        print("槽位、两类整数负载、形状比、长边平方、两种角亏、平均标量曲率",
              n,k4,k6,*[v.str(10) for v in row])

    ctx.prec = 1024
    w,v = 2,1
    previous = None
    for _ in range(40):
        k4,k6 = v,w//2
        assert 3*k4*k4-4*k6*k6 == -1
        if k6>=13:
            n = 1+isqrt(20*k6*k6//3)
            curvature = abs(readings(k4,k6,n)[-1])
            if k6>10**5:
                assert 124<curvature*n**4<126
            if curvature<arb("1e-122"):
                assert previous is not None and previous>arb("1e-122")
                print("Pell 支首次低于所试曲率阈值的槽位与曲率绝对值",n,curvature.str(12))
                break
            previous = curvature
        w,v = 7*w+12*v,4*w+7*v
    else:
        raise AssertionError("所列 Pell 前缀未跨过阈值")


if __name__ == "__main__":
    verify()
