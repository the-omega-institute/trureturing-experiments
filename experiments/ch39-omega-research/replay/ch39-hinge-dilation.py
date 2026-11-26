"""铰链标记数的尺度响应与三阶引力顶点的软腿。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    root = Path(__file__).parent
    flat = run_path(str(root/"ch39-regge-flat-cell.py"))
    local = sum((matrix for _,_,matrix in flat["canonical"]()),sp.zeros(10))
    # 每胞零动量 Hessian 的迹可以直接按二十四个有序单纯形的迹相加。
    assert 24*local.trace() == -7
    s,alpha = sp.symbols("s alpha",positive=True)
    trace = sp.Integer(-7)
    for power,expected in ((0,-2),(1,-1),(alpha,alpha-2)):
        coefficient = s**power*trace/s**2
        assert sp.simplify(sp.diff(coefficient,s).subs(s,1)-expected*trace)==0
    # 在平均数小于半时，比较任何非负整数负载的三阶中心矩与伯努利下界。
    n,mu = sp.symbols("n mu")
    polynomial = (n-mu)**3-((1-n)*(-mu)**3+n*(1-mu)**3)
    assert sp.factor(polynomial) == n*(n-1)*(n+1-3*mu)
    sample = [Fraction(7,10),Fraction(1,5),Fraction(1,10)]
    mean = sum(i*p for i,p in enumerate(sample))
    third = sum((i-mean)**3*p for i,p in enumerate(sample))
    assert mean< Fraction(1,2) and third>mean*(1-mean)*(1-2*mean)
    print("归一零动量 Hessian 的迹、Einstein 与固定标记的软缩放三阶迹",
          trace,-trace,-2*trace)


if __name__ == "__main__":
    verify()
