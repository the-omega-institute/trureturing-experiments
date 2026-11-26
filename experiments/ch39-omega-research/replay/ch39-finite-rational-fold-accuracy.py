"""只用有限次有理运算逼近线性 Higgs 块及其方向导数，同时严格保持全折叠酉性。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import acb_mat, arb, ctx


def bounds(steps):
    q = Fraction(1, 4)
    n = 2**steps
    delta = q**n
    den = 1-q**(n+1)
    dk = n*q**(n-1)/den+(n+1)*q**n/den**2
    dr = 2*delta + (1+q)*dk
    # y(dq/dy)=2q(1+q)/(1−q)，F=yR。所有估计覆盖零背景。
    derivative = delta + 2*q*(1+q)/(1-q)*dr
    logarithmic = 2*q*(1+q)/(1-q)*dr/(1-delta)
    return delta, derivative, logarithmic


def singular_pair(y,steps):
    x = Fraction(1)
    for _ in range(steps):
        x = (x + (1-y*y)/x)/2
    h = y/(1+x)
    return 2*h/(1+h*h), (1-h*h)/(1+h*h)


def real_zero_mixing():
    flavor = run_path(str(Path(__file__).with_name("ch39-chiral-neutral-flavor-phase.py")))
    ctx.prec = 384
    matrix = flavor["reference_matrix"]()
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    scale = arb(4)/(5*arb(".8616"))
    for number in (22027,28325):
        physical = matrix(flavor["singlet"](number)[0])*scale
        assert sum(abs(physical[i,j])**2 for i in range(3) for j in range(3)) < arb(".64")
        square = identity - physical.conjugate().transpose()*physical
        x = identity
        for _ in range(4):
            x = (x + square*x.inv())/2
        h = physical*(identity+x).inv()
        f = 2*h*(identity+h.conjugate().transpose()*h).inv()
        raw, folded = flavor["read"](physical), flavor["read"](f)
        assert all(a.overlaps(b) and (a-b).abs_upper()<arb("1e-26") for a,b in zip(raw,folded))
        print(number, "四步折叠后的 |Vub|、J", folded[2].str(9),folded[3].str(9))


def verify():
    delta, derivative, logarithmic = bounds(4)
    assert delta < Fraction(24,10**11)
    assert derivative < Fraction(21,10**9) and logarithmic < Fraction(21,10**9)
    print("四步的相对质量、绝对导数、对数斜率误差上界",
          *(float(v) for v in (delta,derivative,logarithmic)))
    for i in range(-8,9):
        y = Fraction(i,10)
        f,c = singular_pair(y,4)
        assert f*f+c*c==1 and abs(f-y)<=abs(y)*delta
    # 小尺寸有理矩形例，无需取谱、开方或指定奇异向量。
    T = sp.Matrix([[sp.Rational(1,4),sp.Rational(1,5),0],[0,sp.Rational(1,3),sp.Rational(1,7)]])
    il,ir = sp.eye(T.cols),sp.eye(T.rows)
    A = il-T.T*T
    X = il
    for _ in range(4):
        X = (X + A*X.inv())/2
    H = T*(il+X).inv()
    fl,fr = (il+H.T*H).inv(),(ir+H*H.T).inv()
    F = 2*H*fl
    U = F.row_join((ir-H*H.T)*fr).col_join(((il-H.T*H)*fl).row_join(-F.T))
    assert U.T*U==sp.eye(5) and U*U.T==sp.eye(5)
    assert (F-T).norm()**2 <= sp.Rational(delta.numerator,delta.denominator)**2*T.norm()**2
    # 表 12 的正耦合拟合与最低阶比值进行中心尺度比较。
    ratio = Fraction(291467,100000)
    assert 1/ratio < Fraction(586,1000)**2
    real_zero_mixing()


if __name__=="__main__":
    verify()
