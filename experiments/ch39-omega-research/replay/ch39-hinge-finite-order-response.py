"""正的有限标记规则和 Regge 软顶点的逐阶尺度响应。"""
from fractions import Fraction
from math import comb, factorial

import sympy as sp


def expected_marks(squared_scale, order):
    occupancy = squared_scale/(1+squared_scale)
    contrast = 2*occupancy-1
    return 2*occupancy*sum(Fraction(comb(2*j,j),4**j)*contrast**(2*j)
                           for j in range(order+1))


def verify():
    x = sp.Symbol("x",positive=True)
    for order in range(5):
        occupancy = x/(1+x)
        m = expected_marks(x,order)
        # 将面积平方反射 x -> 1/x 后的恒等式，以及在参考点的归一。
        assert sp.cancel(m-x*m.subs(x,1/x)) == 0
        assert m.subs(x,1) == 1
        assert sp.cancel(sp.diff(m,x).subs(x,1)-sp.Rational(1,2)) == 0
        if order == 0:
            assert sp.cancel(2*occupancy-m) == 0
        else:
            assert m.subs(x,4)>2*occupancy.subs(x,4)
        s = sp.Symbol("s")
        hessian = sp.cancel(m.subs(x,s**2)/s**2)
        degree = 2*order+2
        # 在 s=1 处除去已知阶数的零，精确比较下一项。
        factor = sp.cancel((hessian-1/s)/(s-1)**degree)
        assert factor.subs(s,1) == -sp.Rational(comb(2*order+2,order+1),4**(order+1))
        print("每胞候选数、匹配的最高软腿数、首个差异的导数系数",
              2*(order+1),degree-1,factor.subs(s,1)*factorial(degree))
    # 单次激活的角源概率不超过 1，标记总数有确定的有限上界。
    for ratio in (Fraction(1,9),Fraction(1),Fraction(9),Fraction(100)):
        for order in range(8):
            mean = expected_marks(ratio,order)
            assert 0<mean<=2*(order+1) and mean*mean<=ratio


if __name__ == "__main__":
    verify()
