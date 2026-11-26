"""两个等深电弱真空的最低次势、单态实现与同一场表的端点配平。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


broken = run_path(str(Path(__file__).with_name("ch39-broken-electroweak-vacuum.py")))
graded = broken["graded"]


def schur_example(v, h):
    t = (h/v)**2
    g = arb(2)
    lam = t * (g-t)/(2*(g-2*t))
    other = v * (2*lam*g/t).sqrt()
    slope = (2*lam-t)/(2*lam)
    cosine = 1/(1+slope**2).sqrt()
    sine = slope*cosine
    assert ((2*lam-t)*cosine-2*lam*sine).contains(0)
    assert ((2*lam+g-t)*sine-2*lam*cosine).contains(0)
    triple = 6*v*(lam*(cosine-sine)*cosine**2+g*cosine*sine**2)
    quartic = 24*(lam*cosine**4/4+g*cosine**2*sine**2/2)
    assert arb("377.37") < other < arb("377.38")
    print("单态例子的 λ、额外质量、双重态分量、三次与四次接触比分别为",
          *(a.str(10) for a in (lam,other,cosine,triple*v/(3*h**2),quartic*v**2/(3*h**2))))
    return lam, g, other


def endpoints(light, a, v, h, planck, lam, g, other):
    # 对称点无质量层 24−90，双重态四个实分量等质量，单态另有一质量。
    symmetric = [(-66, arb(0), 0), (4, g*(v/2/planck)**2, 0),
                 (1, 2*lam*(v/planck)**2, 0)]
    massive = [(net,z,-4 if z.is_zero() else 0) for _,net,_,z in light]
    massive += [(1,(other/planck)**2,0),(0,a,12),(0,2*a,-12),(0,3*a,4)]
    def counts(levels):
        assert all(not z.overlaps(w) for j,(_,z,_) in enumerate(levels) for _,w,_ in levels[:j])
        pairs = [broken["minimal_partners"](net,goal) for net,_,goal in levels]
        return tuple(map(sum,zip(*pairs)))
    assert counts(symmetric)==(69,2) and counts(massive)==(109,12)
    # 十组相同质量、每组四标量加一 Dirac，可在对称点用完相同的新增场。
    assert 109-69 == 4*(12-2)
    return broken["light_completion"](light,a)+[(arb(1),(other/planck)**2)]


def shifted_minimum(v, h, planck, M, x):
    # 固定引力截止，v 处以一圈 K=1 标定。求 U=树势+四层热核的破缺局部极小。
    eps = arb("1e-54")
    r = (v*(1-eps)).union(v)
    y = (r**2/(M*planck))**2/x
    i0 = (9*arb(3).log()-12*arb(2).log())/2
    integ = i0-arb(0).union(y.upper())
    scale = 1/(8*arb.pi()**2*M**4)
    # 对 I(y) 的导数用 (1−exp(−y))/y=∫₀¹exp(−ay)da。
    first = -arb(0).union(arb(1))
    second = arb(0).union(arb(3)/2)
    grad = scale*r**7*(8*integ+4*y*first)
    curvature = scale*r**6*(56*integ+76*y*first+16*y**2*second)
    # 树势导数以 t=r/v 表示。端点导数用 −eps(2−eps) 避免相消。
    t = r/v
    tree_hessian = h**2*(30*t**4-24*t**2+2)/8
    t0 = 1-eps
    left_grad = h**2*v*t0*(-eps*(2-eps))*(3*t0**2-1)/4
    assert left_grad+grad < 0 < grad and tree_hessian+curvature > h**2/2
    dx = grad/(tree_hessian+curvature)
    assert dx>0 and dx/v < arb("8e-56")
    # 下降量 ≤ max F′(v)·(v−r*)，因为区间上 U′ 单调。
    rho = scale*v**8*integ
    assert grad*dx/rho < arb("7e-55")
    return dx/v


def verify():
    chiral = graded["minimal"]["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec=512
    source=chiral["source"]["spectral_supplier"](model)
    v, h, planck = arb("246.22"), arb("125.25"), arb("1.22089e19")
    light = broken["light_spectrum"](planck)
    lam,g,other = schur_example(v,h)
    print("微扰单态质量上界 GeV、六次势垒 GeV⁴、树级壁张力 GeV³",
          *((item).str(11) for item in (arb(2).sqrt()*4*arb.pi()*v**2/h,h**2*v**2/54,h*v**2/8)))
    for number in (22027,28325):
        src = source(acb.zeta_zero(number).imag)
        bands = {q:chiral["ordered"]([h for h,charge in src if charge==q])
                 for q in (0,*model["charges"])}
        heavy = graded["minimal"]["fields"](bands)
        paired = graded["paired_spectrum"](heavy)[1:]
        M = 32*arb("2.435e18")/724**2*bands[2][15]
        a = (v**2/(M*planck))**2
        x_six = broken["match"](paired+broken["light_completion"](light,a))
        dx = shifted_minimum(v,h,planck,M,x_six)
        x = broken["match"](paired+endpoints(light,a,v,h,planck,lam,g,other))
        i0 = (9*arb(3).log()-12*arb(2).log())/2
        ratio = a**2*i0/(arb.pi()*arb("1e-122"))*(1-arb(0).union(a/(x*i0)))
        print(number,"六次势的相对驻点偏移、单态模型端点的 Q/M_P 与曲率目标比",
              dx.str(10),x.sqrt().str(10),ratio.str(12))


if __name__=="__main__":
    verify()
