"""重壳对胶子振幅的容量与顶质量层的额外场数。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
higgs = run_path(str(root / "ch39-heat-top-higgs-gluons.py"))
chiral = higgs["broken"]["graded"]["minimal"]["chiral"]


def local_moments():
    # 质量平方对 u=r²/v² 的斜率，以 z_t(v) 为单位。轻重颜色此处只在顶层计。
    # k=0 时均值和方差已强制十二个正向斜率等于顶斜率一。
    # k=1：一个辅助 Dirac 和四个实标量。
    d = (-15-arb(105).sqrt())/6
    m = d+3
    radius = (-6*(d+1)).sqrt()
    slopes = [m-radius]*2 + [m+radius]*2
    assert d < -4 and all(x.is_finite() for x in slopes)
    for j in range(1,4):
        assert (sum(x**j for x in slopes)-12-4*d**j).contains(0)
    assert abs(sum(x**4 for x in slopes)-12-4*d**4) > 1
    # 对线性质量平方，验证全背景构造：十二个固定色标量由三个固定 Dirac 补齐，
    # 真实顶的十二个负偏振由十二个随顶变化的规范单态标量补齐。
    for u in (Fraction(1,2),Fraction(1),Fraction(3,2)):
        net = {}
        for z,n in [(1,12),(u,-12),(1,-4*3),(u,12)]:
            net[z] = net.get(z,0)+n
        assert all(n==0 for n in net.values())
    print("一点三阶热迹匹配：附加 Dirac 质量平方斜率与两组标量斜率",
          d.str(8), slopes[0].str(8), slopes[-1].str(8))


def capacity(number, model, source):
    values = source(acb.zeta_zero(number).imag)
    bands = {q:chiral["ordered"]([h for h,z in values if z==q]) for q in (0,*model["charges"])}
    unit = arb("2.435e18")/724**2
    v,mh = higgs["v"],higgs["mh"]
    masses = []
    for q,index in ((1,2),(2,1),(4,1)):
        assert len(bands[q][3:]) == 16
        masses.extend((unit*(1 if q==1 else 32)*h,index,"f") for h in bands[q][3:])
    assert len(bands[0]) == 18
    masses.extend((unit*32*h,3,"V") for h in bands[0])
    total = arb(0)
    for mass,index,spin in masses:
        y = (mh/(2*mass))**2
        assert 0 < y < arb("1e-15")
        # (asin √y)²=y+y²/3+y³R，0≤R≤1/(1−y)，由正系数的相邻系数比小于一。
        R = arb(0).union(1/(1-y))
        if spin=="f":
            amplitude = 2*(arb(2)/3 + y/3 + (y-1)*y*R)
        else:
            amplitude = -(7 + 2*y + 3*(2*y-1)*y*R)
        total += index*v/mass*abs(amplitude)
    info = higgs["broken"]["light_spectrum"](arb(1))
    sm = sum(higgs["loop"](z.sqrt())[0] for name,_,_,z in info if name in ("u","d","s","c","b","t"))
    needed = abs(sm+5*higgs["loop"](arb("172.5"))[1])-arb("1.10")*abs(sm)
    remaining = needed*(1-(arb("4.1")/arb("38.25")).sqrt())
    assert remaining > 1
    print(number,"重壳每单位最大质量斜率的振幅上界、所需最大质量斜率下界",
          total.upper().str(9),(remaining/total).lower().str(9))
    assert total<arb("1e-8") and remaining/total>10**8


def verify():
    local_moments()
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 224
    source = chiral["source"]["spectral_supplier"](model)
    for number in (22027,28325):
        capacity(number, model, source)


if __name__=="__main__":
    verify()
