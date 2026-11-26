"""在重标量对电弱背景响应的模型中，同时匹配引力、规范交汇和径向 Higgs 驻点。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def verify(root):
    crossing = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    heat = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    ctx.prec = 256
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    quartic = arb(".129")
    positive = arb(".9108").union(arb(".9110"))
    negative = -positive
    portals = (negative, positive)
    shifts = [arb(0)]*14 + [rho*p for p in portals]
    # 整个门户矩形上，每个候选点使用相应移质谱和自己的精确线性系统解。
    sites, weights, rows = crossing["verify"](root, return_measure=True, scalar_shifts=shifts, check_path=True)
    masses = [row[2] for row in rows[49:65]][-2:]
    pi = arb.pi()

    def coefficients(z):
        i = sum(c*x*heat(2, z/x) for c, x in zip(weights, sites, strict=True))
        j = sum(c*heat(1, z/x) for c, x in zip(weights, sites, strict=True))
        return i/(8*pi*pi), -j/(16*pi*pi)

    (m15, b15), (m16, b16) = [coefficients(z+rho*p) for z, p in zip(masses, portals, strict=True)]
    assert m15 < 0 and m16 < 0 and b15 > 0 and b16 > 0
    ratio = m16/m15
    mapped_pos = (quartic/(b16+b15*ratio**2)).sqrt()
    mapped_neg = -ratio*mapped_pos
    assert arb(".9109025") < mapped_pos < arb(".9109027")
    assert -arb(".9109026") < mapped_neg < -arb(".9109024")
    print("整矩形的门户自映射，第 15、16 层", mapped_neg.str(15), mapped_pos.str(15))

    # 固定不动点处的同一组五权，以三阶热核导数包住 W(0)-W(ρv)。
    rho_segment = arb(0).union(rho)
    abs_third = sum(abs(p)**3/(z+rho_segment*p).lower() for z, p in zip(masses, portals, strict=True))
    abs_third *= sum(abs(w) for w in weights)/(8*pi*pi)
    assert abs_third < arb("1e10")
    central = 8*pi*quartic*rho**2
    error = 8*pi*abs_third*rho**3/6
    difference = central+arb(0, error.upper())
    original = arb("1e-122")+difference
    print("对称原点与电弱驻点的 W 差、余项、原点 W 相对目标倍数",
          difference.str(23), error.upper().str(12), (original/arb("1e-122")).str(21))
    assert arb("1.34e-67") < difference < arb("1.35e-67")
    assert error/central < arb("1e-23")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
