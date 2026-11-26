"""固定五原子截止时，背景 Feynman 规范下轻 W、Z 和顶夸克的电弱势曲率。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root):
    crossing = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    sites, weights, rows = crossing["verify"](root, return_measure=True)
    heat = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    prefactor = 1/(8*arb.pi()**2)
    constant = sum(c*x.log() for c, x in zip(weights, sites, strict=True))-arb.const_euler()
    remainder = sum(abs(c)/x for c, x in zip(weights, sites, strict=True))

    def activation(z):
        return sum(c*x*heat(2, z/x) for c, x in zip(weights, sites, strict=True))

    def light(degrees, slope):
        z = slope*rho
        assert z*remainder < arb("1e-19")
        j = constant-z.log()+arb(0, (z*remainder).upper())
        return degrees*prefactor*slope*activation(z)/4, -degrees*prefactor*slope**2*j/4

    def spectrum(g1, g2, yt):
        # g1 为大统一归一超荷；三类偏振含 W+、W- 两个实生成元及三种顶颜色。
        kw = g2**2/2
        kz = (g2**2+3*g1**2/5)/2
        return light(6, kw), light(3, kz), light(-12, yt**2)

    def heavy(z):
        j = sum(c*heat(1, z/x) for c, x in zip(weights, sites, strict=True))
        return prefactor*activation(z), -prefactor*j

    masses = [z for _, _, z in rows[49:65]]
    unit = arb(-1).union(arb(1))
    aa, bb = zip(*(heavy(z+unit*rho) for z in masses), strict=True)
    assert all(a < 0 < b for a, b in zip(aa, bb, strict=True))
    denominator = sum(a**2/b for a, b in zip(aa, bb, strict=True))
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    w, z, top = spectrum(g1, g2, yt)
    assert arb("-.00141574") < w[0] < arb("-.00141571")
    assert arb("-.00093580") < z[0] < arb("-.00093577")
    assert arb("-.03332010") < w[1] < arb("-.03332006")
    assert arb("-.02892379") < z[1] < arb("-.02892375")
    a, b = [sum(axis) for axis in zip(w, z, top, strict=True)]
    assert arb(".00797274") < a < arb(".00797275") and arb(".796672") < b < arb(".796674")
    lower = b+a*a/denominator
    assert lower > arb(".823896")
    print("W、Z、顶的各项一二阶导", *[(x.str(13), y.str(13)) for x, y in (w,z,top)])
    print("合计轻场斜率、曲率、完整驻点的曲率下界", a.str(15), b.str(15), lower.str(15))

    a0, b0 = zip(*(heavy(z0) for z0 in masses), strict=True)
    direction = [-m/k for m, k in zip(a0, b0, strict=True)]
    mid = (a/denominator).mid()
    left, right = mid*(1-arb("1e-6")), mid*(1+arb("1e-6"))

    def total(t):
        slope, curvature = a, b
        for z0, u in zip(masses, direction, strict=True):
            lam = t*u
            m, k = heavy(z0+lam*rho)
            slope += m*lam
            curvature += k*lam**2
        return slope, curvature

    assert total(left)[0] > 0 > total(right)[0]
    domain = left.union(right)
    witness = total(domain)[1]
    assert witness < arb(".823898")
    portals = [domain*u for u in direction]
    assert all(arb(".1001") < p < arb(".1184") for p in portals)
    print("一阶导精确为零的见证门户范围、驻点曲率", min(p.lower() for p in portals).str(12),
          max(p.upper() for p in portals).str(12), witness.str(15))

    # 同一来源的三个边际区间作独立矩形，仅作输入扰动检验。
    pieces = spectrum(g1+arb(0,".000027"), g2+arb(0,".00003"), yt+arb(0,".0043"))
    s, c = [sum(axis) for axis in zip(*pieces, strict=True)]
    margin = c+s*s/denominator
    assert margin > arb(".805")
    print("三项输入边际矩形上的驻点曲率下界", margin.lower().str(15))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
