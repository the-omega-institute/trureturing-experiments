"""用四份相同二态的有理投影检验铰链占据，核对有限精度的构造。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    flat = run_path(str(Path(__file__).with_name("ch39-regge-flat-cell.py")))
    shapes = {}
    for i,j,k in combinations(range(5),3):
        def length(a,b):
            return Fraction((b-a)*(5-b+a),5)
        x,y,z = length(i,j),length(j,k),length(i,k)
        area_squared = (2*(x*y+y*z+z*x)-x*x-y*y-z*z)/16
        shapes[area_squared] = shapes.get(area_squared,0)+1
    assert set(shapes) == {Fraction(1,5),Fraction(3,20)}
    assert flat["R"](3,20)/flat["R"](1,5)==sp.Rational(3,4)
    # 投影仅在 |0000>、|1111> 两个端点有非零矩阵元。
    w = sp.Matrix([3,-4])
    projector = w*w.T/25
    assert projector.T == projector and projector**2 == projector
    p = sp.Symbol("p")
    target = 2*sp.sqrt(3)-3
    assert sp.expand(3*(1-p)**2-4*p*p-(target-p)*(p+target+6)) == 0
    area_sum = 30/sp.sqrt(5)+20*sp.sqrt(sp.Rational(3,20))
    assert sp.simplify(1/area_sum-sp.sqrt(5)/(30+10*sp.sqrt(3))) == 0
    # 若 d=2 mod 4，所需高斯有理数范数在 3 的赋值为奇数。
    assert {(a*a+b*b)%3 for a in range(3) for b in range(3) if a or b} == {1,2}
    z = Fraction(1)
    for _ in range(7):
        z = (3*z + Fraction(3,4)/z**3)/4
    probability = z*z/(1+z*z)
    residual = abs(3-6*probability-probability**2)
    assert residual/6 < Fraction(1,10**122)
    print("有理投影、相同输入份数、七步后六价占据概率的误差上界",
          projector,4,"1e-122")


if __name__ == "__main__":
    verify()
