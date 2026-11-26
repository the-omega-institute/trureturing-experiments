"""有限三标记读数保持前两矩并使零点选择角亏三次作用。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
import sympy as sp


def verify():
    m,z,eta = sp.symbols("m z eta",real=True)
    tau = eta*(m*(1-m))**3
    pgf = sp.expand((1-m+m*z)**3+tau*(z-1)**3)
    base = (1-m+m*z)**3
    for order in range(3):
        assert sp.simplify((pgf-base).diff(z,order).subs(z,1))==0
    assert sp.expand((pgf-base).diff(z,3)-6*tau)==0
    for k in range(4):
        nominal = comb(3,k)*m**k*(1-m)**(3-k)
        ratio = sp.cancel(comb(3,k)*(m*(1-m))**3/nominal)
        assert sp.simplify(ratio-m**(3-k)*(1-m)**k)==0
    a,b = sp.symbols("a b",positive=True)
    seed = (1+a**2-b**2)/2
    assert sp.factor(seed-a-(1-a-b)*(1-a+b)/2)==0
    # 两次迭代的误差含 (1-a-b)^4；只需验证一次更新的恒等式。
    assert sp.factor(a-2*m*a*a/(m*m+a*a)-a*(m-a)**2/(m*m+a*a))==0
    shell=run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    ctx.prec=192
    eta_values=[shell["phases"](acb.zeta_zero(n).imag)[724].real for n in (11366,11431)]
    b=2*arb(3).sqrt()-3
    a=1-b
    d4=arb(18)/(5*8*arb.pi())
    d6=-arb(3).sqrt()*d4
    difference=-2*(eta_values[0]-eta_values[1])*(a*b)**3*(30*d4**3+20*d6**3)
    assert difference>0
    print("互补平直胞形状变分的曲率三次项，11366 减 11431",difference.str(14))


if __name__=="__main__":
    verify()
