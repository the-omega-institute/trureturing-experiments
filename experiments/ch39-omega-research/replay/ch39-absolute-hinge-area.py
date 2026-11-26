"""用相同副本的计数读数选择绝对面积占据，并检验最低次数。"""
from math import comb

import sympy as sp


def bernstein_mean(weights, probability):
    copies = len(weights)-1
    return sp.expand(sum(value*comb(copies,k)*probability**k*(1-probability)**(copies-k)
                         for k,value in enumerate(weights)))


def verify():
    p,x,r = sp.symbols("p x r", real=True)
    three = [sp.Rational(2,3)+(int(k==3)-x*k)/3 for k in range(4)]
    four = [(1+4*(int(k==4)-x*comb(k,2)/3+x**2))/5 for k in range(5)]
    cubic = bernstein_mean(three,p)
    quartic = bernstein_mean(four,p)
    assert sp.expand(cubic-(sp.Rational(2,3)+p**3/3-x*p)) == 0
    assert sp.expand(quartic-(1+4*(p**2-x)**2)/5) == 0
    for weight in three+four:
        for boundary in (weight,1-weight):
            # 二次多项式在闭区间的端点或临界点取极小。
            points = [0,1]+[v for v in sp.solve(sp.diff(boundary,x),x) if 0<v<1]
            assert all(boundary.subs(x,v)>=0 for v in points)
    lowest = sp.expand(cubic.subs({p:r,x:r*r}))
    assert sp.expand(cubic.subs(x,r*r)-lowest-(p-r)**2*(p+2*r)/3) == 0
    assert sp.expand(5*quartic.subs(x,r*r)-1-4*(p-r)**2*(p+r)**2) == 0
    coefficients = sp.QQ.frac_field(x)
    root_polynomial = sp.Poly(p**2-x,p,domain=coefficients)
    assert root_polynomial.is_irreducible
    assert sp.diff(cubic,p) == root_polynomial.as_expr()
    assert sp.expand(quartic-sp.Rational(1,5)-sp.Rational(4,5)*root_polynomial.as_expr()**2)==0
    print("绝对面积占据的副本下界、具有有理最低代价时的副本下界",3,4)
    print("面积比 r 处的最小读数",lowest,sp.Rational(1,5))


if __name__ == "__main__":
    verify()
