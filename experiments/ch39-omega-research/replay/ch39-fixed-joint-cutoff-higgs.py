"""固定全场五截止解时，重 Higgs 双重态对邻域势的两阶阈值。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root):
    crossing = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    sites, weights, rows = crossing["verify"](root, return_measure=True)
    heat = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    pi = arb.pi()
    masses = [z for _, _, z in rows[49:65]]
    assert len(masses) == 16 and min(z.lower() for z in masses) > arb("1.61822e-11")

    def coefficients(z):
        # V/M_P^4 的 ρ 和 ρ² 系数。每个复标量双重态有四个实偏振。
        i = sum(c*x*heat(2, z/x) for c, x in zip(weights, sites, strict=True))
        j = sum(c*heat(1, z/x) for c, x in zip(weights, sites, strict=True))
        return i/(8*pi*pi), -j/(16*pi*pi)

    linear, quadratic = zip(*(coefficients(z) for z in masses), strict=True)
    assert all(m < 0 < b for m, b in zip(linear, quadratic, strict=True))
    m, b = sum(linear), sum(quadratic)
    w1, w2 = 8*pi*m, 16*pi*b
    assert arb("-1.864172") < w1 < arb("-1.864170")
    assert arb("59.32082") < w2 < arb("59.32086")
    assert arb("-.07417302") < m < arb("-.07417300")
    assert arb("1.1801505") < b < arb("1.1801509")

    quartic = arb(".129")
    v = arb("246.22")/arb("1.22089e19")
    delta = quartic*v*v
    common_mass = -delta/m
    common_quartic = (quartic/b).sqrt()
    assert arb("7.07354e-34") < common_mass < arb("7.07356e-34")
    assert arb("5.9048e-67") < b*common_mass**2 < arb("5.9050e-67")
    assert arb(".330617") < common_quartic < arb(".330619")
    assert abs(m)*common_quartic-delta > arb(".02452290")
    print("单位共同门户 W 的一二阶导、物理质量平方与四次系数",
          *[x.str(17) for x in (w1, w2, m, b)])
    print("质量拟合时共同门户与四次系数、四次拟合所需门户绝对值与质量偏移",
          *[x.str(13) for x in (common_mass, b*common_mass**2, common_quartic, abs(m)*common_quartic)])

    # 奇偶两组抵消线性项，再加精确的极小共同偏移，使净量等于 -δ。
    ratio = sum(linear[::2])/sum(linear[1::2])
    direction = [arb(1) if j%2 == 0 else -ratio for j in range(16)]
    b1 = sum(y*u for y, u in zip(quadratic, direction, strict=True))
    b2 = sum(y*u*u for y, u in zip(quadratic, direction, strict=True))
    offset = common_mass
    size = (-offset*b1 + (offset**2*b1**2 + b2*(quartic-offset**2*b)).sqrt())/b2
    odd, even = offset+size, offset-size*ratio
    assert 0 < odd < arb(".330619") and -arb(".330619") < even < 0
    print("奇偶两组的精确公式所得门户", odd.str(14), even.str(14))

    # 对给定线性项绝对和，两种符号的四次响应分别由各自最大效率行达到。
    efficiency = [y/(x*x) for x, y in zip(linear, quadratic, strict=True)]
    first, second = efficiency[-1], efficiency[-2]
    assert first > second > max(e.upper() for e in efficiency[:-2])
    total = ((2*(quartic*(first+second)-delta**2*first*second).sqrt()
              -delta*(first-second))/(first+second))
    l16 = (total+delta)/(-2*linear[-1])
    l15 = (total-delta)/(2*linear[-2])
    assert arb(".00844554") < total < arb(".00844556")
    assert arb("6.21233e-33") < delta/total < arb("6.21236e-33")
    assert 0 < l16 < 1 and -1 < l15 < 0
    print("最小绝对质量贡献、净量占比及达到它的第15、16层门户",
          *[x.str(15) for x in (total, delta/total, l15, l16)])

    # 在实际电弱背景直接匹配一阶驻点与二阶径向曲率。用 Brouwer 的闭矩形自映射。
    rho = v*v/2
    unit = arb(-1).union(arb(1))
    assert all(coefficients(z+rho*unit)[0] < 0 for z in masses)
    positive = arb(".9108").union(arb(".9110"))
    negative = -positive
    m15, b15 = coefficients(masses[-2]+rho*negative)
    m16, b16 = coefficients(masses[-1]+rho*positive)
    ratio = m16/m15
    mapped_positive = (quartic/(b16+b15*ratio**2)).sqrt()
    mapped_negative = -ratio*mapped_positive
    assert mapped_positive > arb(".9108") and mapped_positive < arb(".9110")
    assert mapped_negative > -arb(".9110") and mapped_negative < -arb(".9108")
    print("实际电弱驻点自映射的输出包围", mapped_negative.str(14), mapped_positive.str(14))

    # 各质量在 [0,rho] 维持正值。热核三阶导的界也适用于任意裸二次、四次势。
    variation = sum(abs(c) for c in weights)
    assert variation < 4
    def third_bound(limit):
        mass_floor = min(z.lower() for z in masses) - limit*rho
        assert mass_floor > 0
        return len(masses)*limit**3*variation/(8*pi**2*mass_floor)

    # 从驻点向对称原点作二阶 Taylor 展开。
    error = third_bound(arb(1))*rho**3/6
    vacuum_drop = 8*pi*(quartic*rho**2 + arb(0, error.upper()))
    assert arb("1.34077e-67") < vacuum_drop < arb("1.34079e-67")
    assert third_bound(arb(40_000_000)) < 6*quartic/rho
    print("电弱驻点降低的 W、两处 W 相等所需三阶导的绝对值",
          vacuum_drop.str(14), (6*quartic/rho).str(14))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
