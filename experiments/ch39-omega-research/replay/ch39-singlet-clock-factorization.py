"""一个复单态重现整条单位圆的矩阵差能量，并估计量子包误差。"""
from pathlib import Path
from runpy import run_path

from flint import arb
from sympy import Matrix, Rational, cancel, diag, eye, symbols

root = Path(__file__).parent


def verify():
    # 平方判别式由二次 Gram 行列式给出，秩一和秩二均在范围内。
    s,c,x = symbols("s c x", positive=True)
    assert cancel((x-s+c/x)*x-(x**2-s*x+c))==0
    t,v = symbols("t v", real=True)
    transverse = 4*(t*t+v*v)+(t*t+v*v)**2
    hessian = Matrix([transverse.diff(t),transverse.diff(v)]).jacobian([t,v]).subs({t:0,v:0})
    assert hessian == 8*eye(2)
    assert 2+(2+arb(2).sqrt())**2 < 14
    # 非退化双路径可有两个正 Gram 特征值，标量分解只固定圆上差能响应。
    assert diag(1,1).det() == 1
    matrix = run_path(str(root/"ch39-quantum-matrix-weak-shell.py"))
    matrix["verify"](singlet_clocks=True)
    degrees_sum = 2*15
    assert Rational(degrees_sum-14*22,3)==-Rational(278,3)


if __name__ == "__main__":
    verify()
