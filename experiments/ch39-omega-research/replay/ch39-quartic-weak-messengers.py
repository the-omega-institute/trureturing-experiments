"""三个双基本中介场给出最高四次的标量势及其质量、耦合代价。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import cancel, symbols

root = Path(__file__).parent
q = arb


def verify():
    l,mu,r,s = symbols("l mu r s", positive=True)
    density = 2*l*(r-1)**2+mu*(r-s)**2
    trial = (2*l+mu*s)/(2*l+mu)
    assert cancel(density.subs(r,trial)-2*l*mu*(s-1)**2/(2*l+mu))==0
    ctx.prec = 128
    p = 23226
    heat = q("5e-18")
    source = run_path(str(root/"ch39-thermal-plaquette-occupation.py"))
    joint = source["joint"]
    theta = joint["winding"]["angle"](11431)
    alpha = 2*q.pi()/3
    delta = alpha-p*theta
    chord = lambda u: 2*(1-u.cos())
    e = chord(delta)/p
    # mu=8, 4<=l<=6。两个额外的平方配方以串联方式给下界。
    inverse = q(1)/8+q(1)/8
    slope = 3/(1+3*inverse)**2
    correction = 3*inverse*e**2
    early = slope*4*(theta/2).sin()*min((delta+theta/2).sin().lower(),alpha.sin().lower())/p-q("2e-16")
    middle = slope*(chord(theta-delta)/(p+21000)-e)
    late = 21000*q("1e-17")-3*e
    count_gap = min(early.lower()-correction,middle.lower()-correction,late.lower())/heat-q(21010).log()
    assert count_gap>12000
    # 两枢纽的站点度由三个中介各增加一，其余各站点亦有负一圈系数。
    assert 3+3-22==-16 and 2*(15+3)-14*22==-272
    a_low = 1/(1+e/8)
    s_low = 1/(1+e/4)
    tolerance = q("1e-9")
    inertia_upper = ((1+tolerance)**2-1)/(2*s_low**2)
    assert inertia_upper < q("1.001e-9")
    assert ((1+2*q("5e-10"))/a_low).sqrt()-1 < tolerance
    assert 4*q.pi()*inertia_upper**2 < q("1.26e-17")
    print("经典矩阵极小、量子链内振子下错长度指数下界",count_gap.str(12),
          "1 ppb 对中介动能比的必要上界",inertia_upper.str(12),
          "树级标量接触条件对 g 的必要上界",(4*q.pi()*inertia_upper**2).str(12))


if __name__ == "__main__":
    verify()
