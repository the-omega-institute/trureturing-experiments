"""Finite nonminimal-coupling shifts and the scalar-exchange curvature coefficient."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Rational, diff, expand, fraction, symbols


def verify():
    ctx.prec = 160
    root = Path(__file__).parent
    graph = run_path(str(root/"ch39-weighted-arc-zero-phase.py"))
    a, _, fit = run_path(str(root/"ch39-shell-zero-mass-hierarchy.py"))["parameters"]()
    num, den = fraction(graph["response"](3))
    target = a("1.2209e19")
    for j in (7002,7003,7004,11485,11490):
        h = acb.zeta_zero(j).imag
        u = 1-1/(2*(h**2+a("1/4")))
        e = graph["evaluated"](num,u)/graph["evaluated"](den,u)
        low,v,_,_,_ = fit(e,a(1))
        m2 = (v/low)**2*(2*arb.pi()/(3*arb(5).sqrt()/2))*(h*h+a("1/4")+30+a("2/15"))
        scale = 160*arb.pi()*v*v
        xi = (target**2-m2)/scale
        mass = (m2+scale*a(-1).union(a(1))).sqrt()
        if j == 7002:
            assert mass < a("1.22085e19") and xi > 68
        if j == 7003:
            assert a("1.22085e19") < mass < target and a("40.3143") < xi < a("40.3145")
            assert 944 < (10*xi)**2/172 < 946
        if j == 7004:
            assert mass > a("1.22095e19") and xi < -76
        if j >= 11485:
            assert -692000 < xi < -690000
        print(j,"common nonminimal coupling for exact matching",xi.str(11))
    x,y = symbols("x y",real=True)
    A,B = (1+x)**2,(1+y)**2
    pot = 3*(A-1)**2+7*(B-1)**2+9*(A-1)**2+24*A*(A-1)**2
    assert [diff(pot,x,2).subs(x,0)/2,diff(pot,y,2).subs(y,0)/2,diff(pot,x,y)] == [144,28,0]
    q = x*x/16+7*y*y/4
    assert expand(172*q-(3*x+7*y)**2-Rational(7,4)*(x-12*y)**2) == 0


if __name__ == "__main__":
    verify()
