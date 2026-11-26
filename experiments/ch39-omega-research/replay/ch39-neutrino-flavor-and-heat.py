"""壳质量与 Majorana 味矩阵给出的带电单态衰变，以及破缺真空的补场数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


broken = run_path(str(Path(__file__).with_name("ch39-broken-electroweak-vacuum.py")))
graded = broken["graded"]


def decay_bounds(mu_ev):
    # 正常排序、m1=0。只用有限区间、PMNS 幺正性与三角不等式。
    m2 = arb(".000068").union(arb(".000083")).sqrt()
    m3 = arb(".00238").union(arb(".00265")).sqrt()
    # sin² θ23 在 [0.40,0.63] 上，x(1-x) 的最小值在 0.63 处。
    mu_tau = m3.lower()*arb(".975")*(arb(".63")*arb(".37")).sqrt() - m2.upper()/2
    em = m3.upper()*(arb(".025")*arb(".63")).sqrt() + m2.upper()*arb(".36").sqrt()
    v = arb("246.22")
    me, mm, mtau = map(arb, (".000511",".10566","1.77686"))
    strength = 4*(mtau**2-mm**2)*mu_tau/(v**2*mu_ev)
    eps = ((mm**2-me**2)/(mtau**2-mm**2)*em/mu_tau)**2
    # x_eμ / Σx_ij ≤ eps。最多两个轻子相空间损失的单纯形范围。
    w_min = (1-(mtau/arb("80.38"))**2)**2
    # 只有 eμ 对额外给出一份轻轻子分宽；eτ 或 μτ 对的重通道至少给 w_min。
    # 对所有系数，B_e+B_μ ≤ (1+eps)/(1+eps+w_min*(1-eps))。
    upper = (1+eps)/(1+eps+w_min*(1-eps))
    width = arb("80.38")/(16*arb.pi()) * 2*w_min*strength**2
    life = arb("6.582119569e-25")/width
    assert mu_tau > arb(".018") and eps < arb(".000006")
    assert upper < arb(".50026") and life < arb("5e-19")
    print("μτ 质量元下界/eV，实际轻子相空间的 B_e+B_μ 上界，寿命上界/s",
          mu_tau.str(9), upper.str(9), life.str(9))
    return m2, m3


def rebalanced_heat(bands, mu, spectrum, masses_ev):
    planck = arb("1.22089e19")
    a = (mu/planck)**2
    levels = []
    for name,g,alpha,z in spectrum:
        if z.is_zero():
            g += 2*len(masses_ev)
            alpha -= arb(len(masses_ev))/6
        levels.append((g,alpha,z,-4 if z.is_zero() else 0))
    levels.extend((-2, arb(1)/6, (m*arb("1e-9")/planck)**2, 0) for m in masses_ev)
    levels.extend((0,arb(0),j*a,target) for j,target in ((1,12),(2,-12),(3,4)))
    assert all(not z.overlaps(w) for i,(_,_,z,_) in enumerate(levels)
               for _,_,w,_ in levels[:i])
    rows = graded["paired_spectrum"](graded["minimal"]["fields"](bands))[1:]
    nscalar, ndirac = 896, 162
    for g,alpha,z,target in levels:
        s,d = broken["minimal_partners"](g,target)
        assert g+s-4*d == target
        rows.append((alpha+arb(s)/6+arb(d)/3,z))
        nscalar += s
        ndirac += d
    assert (nscalar,ndirac) == ((1006,174) if len(masses_ev)==2 else (1010,175))
    x = broken["match"](rows)
    # 正四层读数的固有时余项界，与原经热核验证的四层公式相同。
    i0 = (9*arb(3).log()-12*arb(2).log())/2
    ratio = a**2*i0/arb.pi()/arb("1e-122") * (1-arb(0).union(a/(x*i0)))
    print("有质量中微子数、总实标量、总 Dirac、Q/M_P、四层读数/目标",
          len(masses_ev), nscalar, ndirac, x.sqrt().str(9), ratio.str(10))
    return ratio


def verify():
    chiral = graded["minimal"]["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 768
    source = chiral["source"]["spectral_supplier"](model)
    light = broken["light_spectrum"](arb("1.22089e19"))
    for number in (22027,28325):
        src = source(acb.zeta_zero(number).imag)
        bands = {q:chiral["ordered"]([h for h,z in src if z==q]) for q in (0,*model["charges"])}
        heavy = 32*arb("2.435e18")/724**2*bands[2][15]
        mu = arb("246.22")**2/heavy
        print("正零点", number)
        m2,m3 = decay_bounds(mu*10**9)
        r = rebalanced_heat(bands,mu,light,[m2,m3])
        # 另取 m1=0.01 eV，平方质量差沿用同一输入窗口。
        m1 = arb(".01")
        r3 = rebalanced_heat(bands,mu,light,[m1,(m1**2+m2**2).sqrt(),(m1**2+m3**2).sqrt()])
        assert abs(r-1) < arb(".0006") and abs(r3-1) < arb(".0006")


if __name__ == "__main__":
    verify()
