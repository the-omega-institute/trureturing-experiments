"""七式交点附近，带电与中性重标量分量的门户劈裂对弱三分量残差的影响。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq


def slope_enclosure(z, sites, weights, p0, p1, panels):
    # 归一化标量角权；对同一角度，先合并有符号截止，再积质量导数。
    result = arb(0)
    for j in range(panels):
        left = arb(j)/(2*panels)
        right = arb(j+1)/(2*panels)
        h = (left*(1-left)).union(right*(1-right))
        angle = 3*((1-2*right)**2).union((1-2*left)**2)/panels
        a, delta = z+p0*h, (p1-p0)*h
        b = z+p1*h
        exposure = arb(0)
        for x, c in zip(sites, weights, strict=True):
            # exp(-a/x)/a - exp(-b/x)/b。用 expm1 保留端点 h=0 的消去。
            exposure += c*(-a/x).exp()*(delta/(a*b)-(-delta/x).expm1()/b)
        result += angle*exposure
    return result


def verify(root):
    ctx.prec = 192
    t, zero = fmpq(1,2), fmpq(0)
    assert t*(t-1)/3 == -fmpq(1,12) and (-t)*(-t-zero)/3 == fmpq(1,12)
    energy = arb("1.1763969e14").union(arb("1.17639745e14"))
    portal = arb(".42995401871").union(arb(".42995401873"))
    broken = run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"]
    system = broken(root,energy,portal)
    weights = system["weights_at"](portal)[0]
    data = system["original"]
    sites, rho, masses, signs = data["sites"], system["rho"], system["masses"], system["signs"]
    tv = sum(abs(c) for c in weights)
    assert len(masses) == 16 and tv < arb("24.465")
    # B、C 七式根的权重在此盒内。更大的归一化五权族可用统一总变差界。
    lower = [z-rho for z in masses]
    def upper_bound(lam):
        assert min(lower) > lam*rho
        return 25*sum((2*lam*rho/(z-lam*rho)).log1p() for z in lower)/(48*arb.pi())

    unit = upper_bound(arb(1))
    threshold = arb("6.48e21")
    broad = upper_bound(threshold)
    print("七式权总变差、单位门户劈裂界、阈值处劈裂界",tv.str(13),unit.str(15),broad.str(15),flush=True)
    assert unit < arb("1.143e-23") and broad < arb(".07425")
    light = run_path(str(root/"ch39-neutral-residual-all-cutoffs.py"))
    light["verify"]()
    # 此处用该证书已经核验的低端包围及 T <= 25 时的全能标窗余量。
    residual = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    residual += arb(0,arb("4.907e-20"))
    assert residual+arb(0,unit.upper()) < arb("-.07425")
    assert residual+arb(0,broad.upper()) < 0

    max_shift = arb(0)
    direction = []
    for z, s in zip(masses,signs,strict=True):
        zpath = z+s*portal*rho+arb(0,rho.upper())
        derivative = slope_enclosure(zpath,sites,weights,data["p0"],data["p1"],512)
        assert not derivative.contains(0), derivative
        sign = 1 if derivative > 0 else -1
        direction.append(sign)
        max_shift += 2*rho*abs(derivative)/(48*arb.pi())
    assert direction == [1]*16
    assert arb("2.8267e-25") < max_shift < arb("3.0945e-25")
    print("按质量升序，最优中性分量的附加门户号",direction,flush=True)
    print("单位独立分量劈裂的最大正位移下界及上界",
          max_shift.lower().str(17),max_shift.upper().str(17),flush=True)

if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
