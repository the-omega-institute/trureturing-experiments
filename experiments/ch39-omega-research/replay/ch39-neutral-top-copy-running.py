"""顶层中性复制场的门户跑动与逐谱补全后的牛顿标定。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
portal = run_path(str(root/"ch39-portal-yukawa-stability.py"))
broken = run_path(str(root/"ch39-broken-electroweak-vacuum.py"))
graded = broken["graded"]


def running():
    rel = Fraction(1,1000)
    g1, _, g3 = portal["sm"]["hull_inputs"](".467453 .63811 1.0583", rel)
    ys = portal["sm"]["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101", rel)
    hyper = 3*g1**2/5
    gauge = [16*g3**2 + 4*hyper/3]*3 + [16*g3**2 - 2*hyper/3]*3 + [6*hyper]*3
    for i in (0,1,3,4,5):
        # 给全部非顶夸克的每个伙伴最有利的色 Casimir 及超荷扣除。
        gauge[i] -= 20*g3**2 + 2*hyper/3
    p = [y*y for y in ys]
    w = [n*u*u for n,u in zip([12]*6+[4]*3,p)]
    moment = sum(w)
    lower = sum(weight*(u+h) for weight,u,h in zip(w,p,gauge))/moment
    velocity = lower/(8*arb.pi()**2)
    budget = 8*arb.pi()
    eta = lower/budget
    extra = moment*eta**2/2
    print("零额外门户时的斜率下界、8π 四次张量预算下的门户范数比与额外 16π²βλ 下界",
          velocity.lower().str(8), eta.lower().str(8), extra.lower().str(8))
    assert velocity > arb(".235") and eta > arb(".7") and extra > arb("1.6")


def gravity():
    chiral = graded["minimal"]["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 768
    source = chiral["source"]["spectral_supplier"](model)
    planck = arb("1.22089e19")
    spectrum = broken["light_spectrum"](planck)
    ztop = next(z for name,_,_,z in spectrum if name=="t")
    for number in (22027,28325):
        src = source(acb.zeta_zero(number).imag)
        bands = {q:chiral["ordered"]([h for h,z in src if z==q]) for q in (0,*model["charges"])}
        M = 32*arb("2.435e18")/724**2*bands[2][15]
        a = (arb("246.22")**2/(M*planck))**2
        rows = graded["paired_spectrum"](graded["minimal"]["fields"](bands))[1:]
        rows += broken["light_completion"](spectrum,a)
        # 三组中性四标量加一 Dirac，各组给一单位牛顿权，真空净偏振为零。
        x = broken["match"](rows + [(arb(3),ztop)])
        i0 = (9*arb(3).log()-12*arb(2).log())/2
        ratio = a**2*i0/arb.pi()/arb("1e-122")*(1-arb(0).union(a/(x*i0)))
        assert abs(ratio-1)<arb(".0006")
        print(number, "三组额外顶层后的 Q/M_P 和四层真空比值", x.sqrt().str(9), ratio.str(10))


if __name__=="__main__":
    ctx.prec = 224
    running()
    gravity()
