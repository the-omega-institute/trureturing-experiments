"""破缺点真空热迹的二阶抵消对物理 Higgs 衰变到轻标量的约束。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


vacua = run_path(str(Path(__file__).with_name("ch39-two-electroweak-vacua.py")))
broken = vacua["broken"]


def verify_jet_extension():
    # 标准化背景 t=(H†H)/s₀；前三项为所需质量的局部 Taylor 系数，末项选定原点值。
    # 先以 m(t)=sqrt(t) 为例把 m(0) 指定为零。
    c = [Fraction(n, 8) for n in (0, 15, -10, 3)]
    assert sum(c)==1 and sum(i*v for i,v in enumerate(c))==Fraction(1,2)
    assert sum(i*(i-1)*v for i,v in enumerate(c))==Fraction(-1,4)
    assert sum((Fraction(1),Fraction(-1,2),Fraction(-1,8)))==Fraction(3,8)
    # E₂(tz)/t 的导数为 −exp(−tz)/t²；(E₁(tz))′=−exp(−tz)/t。
    for z in (arb(1)/17,arb(3)/2):
        t = arb(5)/4
        e = (-t*z).exp()
        derivative = -z*(t*z).expint(1)/t - (t*z).expint(2)/t**2
        assert (derivative + e/t**2).contains(0)


def width(light, planck, v, mh, bottom):
    widths = {}
    for name, g, _, z in light:
        if g >= 0:
            continue
        m = arb(bottom) if name == "b" else planck*z.sqrt()
        if 2*m >= mh:
            continue
        phase = (1-(2*m/mh)**2).sqrt()
        widths[name] = (-g)*m**4*phase/(8*arb.pi()*v**2*mh)
        sm = (-g)*m**2*mh*phase**3/(32*arb.pi()*v**2)
        assert (widths[name] - sm*4*m**2/(mh**2-4*m**2)).contains(0)
    return widths


def verify():
    ctx.prec=192
    verify_jet_extension()
    planck, v, mh = map(arb,("1.22089e19","246.22","125.25"))
    light = broken["light_spectrum"](planck)
    t = (mh/v)**2
    g = arb(2)
    lam = t*(g-t)/(2*(g-2*t))
    kappa_sq = 1/(1+((2*lam-t)/(2*lam))**2)
    for mass, lo, hi in (("4.18","19.52","19.54"),("2.7","3.70","3.72")):
        widths = width(light,planck,v,mh,mass)
        total = 10**6*sum(widths.values())
        assert arb(lo)<total<arb(hi)
        print("底质量/GeV、底伙伴分宽/keV、全部下界/keV、混合轻态下界/keV",
              mass, (10**6*widths["b"]).str(10), total.str(10), (total*kappa_sq).str(10))


if __name__=="__main__":
    verify()
