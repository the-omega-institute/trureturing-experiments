"""非负有效固有时截止轮廓的曲率障碍和联合编号处的最小峰值。"""
from math import factorial
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


NUMBER = 2_306_888_100_085


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    old = load("chiral-induced-gravity-vacuum")
    lib = load("chiral-minimal-neutral-gravity")
    spectrum = old["spectrum"]
    original, *_ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    broad = spectrum["shell_bands"](original)
    assert old["prefix_cap"](broad, 32) == 0
    chiral = lib["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 256
    raw = chiral["source"]["spectral_supplier"](model)(acb.zeta_zero(NUMBER).imag)
    bands = {q:chiral["ordered"]([h for h,k in raw if k == q]) for q in (0,*model["charges"])}
    rows = lib["fields"](bands)
    assert len(rows) == 99 and rows[0][0] == -62 and rows[0][2].is_zero()
    assert (rows[0][1]-arb(1)/6).contains(0)
    mass_floor = min(z.lower() for _,_,z in rows[1:])
    assert mass_floor > 0
    # 仅在最后的微小相对误差比较中需要高精度。质量带延用已核验的区间。
    ctx.prec = 1408
    pi = arb.pi()
    ceiling = arb(545)/6
    assert (sum(a for _,a,_ in rows if a>0)-ceiling).contains(0)
    curvature = 31*pi/ceiling**2
    assert curvature > arb(".01180")
    k, w = 1/(12*pi), 31/(4*pi)
    epsilon = arb(10)**-122
    cutoff = arb(10)**-60
    u_floor = mass_floor/cutoff
    alpha = sum(abs(a) for _,a,_ in rows[1:])
    spin = sum(abs(g) for g,_,_ in rows[1:])
    heat_error = max((6*alpha).upper(), (arb(spin)/62).upper())*(factorial(8)+factorial(9))/u_floor**8
    error = arb(10)**-330
    assert heat_error < error
    assert error*(1+error) < (1-error)**2
    factor = (1-error).union(1+error)
    unrestricted = 1116*pi/epsilon
    light_t = epsilon/93
    assert (w/k**2/epsilon-unrestricted).contains(0)
    # W/K 在红外严格递增。以 300 位十进相对精度圈住 W(t)=εK(t) 的唯一小正根。
    near = arb(10)**-302
    left, right = light_t*(1-near), light_t*(1+near)
    assert 0 < left < right < cutoff
    assert ((w*left/k)*factor/factor) < epsilon < ((w*right/k)*factor/factor)
    root = left.union(right)
    infimum_ratio = 1/(k*root*factor*unrestricted)
    assert abs(infimum_ratio-1) < arb(10)**-300
    # 对偶切线在红外外始终给负导数。W′≥2wx，K′≤ceiling/(2π)。
    slope = k*factor/(2*w*root*factor)
    assert (ceiling/(2*pi))/(2*w*cutoff) < slope
    assert 1-epsilon*slope > 0

    # 严格正的两截止原子。给 t 选相反符号端点，连续性给精确根。
    r = arb(10)**-25
    center = light_t/(1+r)
    width = arb(10)**-90
    lo, hi = center*(1-width), center*(1+width)
    def vacuum_after_newton(t):
        s = r*t
        delta_w_over_delta_k = (w/k)*(t+s)*factor/factor
        return w*s*s*factor+(1-k*s*factor)*delta_w_over_delta_k
    assert vacuum_after_newton(lo) < epsilon < vacuum_after_newton(hi)
    t = lo.union(hi)
    s = r*t
    peak = (1-k*s*factor)/(k*(t-s)*factor)
    assert 0 < s < t < 1 and peak > 1
    shift = peak/unrestricted-1
    assert 2*r < shift < 2*r+3*r*r
    print("全相位单位高度曲率下界、实际编号、最小重质量平方比下界",
          curvature.str(15), NUMBER, mass_floor.str(15))
    print("红外四项相对误差共同上界、理想峰值常数", heat_error.str(12), (1116*pi).str(17))
    print("最小峰值下确界的相对误差上界、最优端点 x", abs(infimum_ratio-1).upper().str(5), root.str(19))
    print("严格正两截止 s、t、峰值相对理想值的超额", s.str(19), t.str(19), shift.str(29))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
