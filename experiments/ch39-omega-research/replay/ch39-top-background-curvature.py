"""既定五原子截止下三色顶夸克和十六层重标量的 Higgs 电弱曲率。"""
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

    def mass_response(z):
        return sum(c*x*heat(2, z/x) for c, x in zip(weights, sites, strict=True))

    def heavy_coefficients(z):
        j = sum(c*heat(1, z/x) for c, x in zip(weights, sites, strict=True))
        return prefactor*mass_response(z), -prefactor*j

    def top_response(y):
        z = y*y*rho
        # 归一化给出 -log(z) 的精确单位系数；Ein 的余量不超过 z/x。
        assert z*remainder < arb("1e-19")
        j = constant-z.log()+arb(0, (z*remainder).upper())
        t = -3*prefactor*y**2*mass_response(z)
        c = 3*prefactor*y**4*j
        return t, c, j

    masses = [z for _, _, z in rows[49:65]]
    box = arb(-1).union(arb(1))
    a, b = zip(*(heavy_coefficients(z+box*rho) for z in masses), strict=True)
    assert all(x < 0 < w for x, w in zip(a, b, strict=True))
    assert min(w.lower() for w in b) > arb(".1337")
    d = sum(x*x/w for x, w in zip(a, b, strict=True))

    slope, curvature, j = top_response(arb(".8616"))
    floor = curvature+slope*slope/d
    assert arb(".01032425") < slope < arb(".01032426")
    assert arb(".858915") < curvature < arb(".858919")
    assert arb("41.0200") < j < arb("41.0202")
    assert floor > arb(".904566")
    print("1 TeV 中心值的顶一阶导、二阶导、对数核与驻点曲率下界",
          *[x.str(15) for x in (slope, curvature, j, floor)])

    # 给出几乎达到驻点下界的确切门户存在证书。
    a0, b0 = zip(*(heavy_coefficients(z) for z in masses), strict=True)
    direction = [-x/w for x, w in zip(a0, b0, strict=True)]
    center = (slope/d).mid()
    radius = center*arb("1e-6")
    left, right = center-radius, center+radius
    def coupled(t):
        total_slope, total_curvature = slope, curvature
        for z, u in zip(masses, direction, strict=True):
            portal = t*u
            x, w = heavy_coefficients(z+portal*rho)
            total_slope += portal*x
            total_curvature += portal**2*w
        return total_slope, total_curvature
    assert coupled(left)[0] > 0 > coupled(right)[0]
    portals = [left.union(right)*u for u in direction]
    assert all(0 < p < arb(".16") for p in portals)
    witness = coupled(left.union(right))[1]
    assert witness < arb(".904570")
    print("一阶精确消去见证的十六层门户范围、二阶导", min(p.lower() for p in portals).str(10),
          max(p.upper() for p in portals).str(10), witness.str(15))

    # 表列边际输入整窗，顶部对数单独处理以保留跨五原子的抵消。
    ends = [arb(".8573") + arb(".0086")*i/8 for i in range(9)]
    pieces = [top_response(lo.union(hi)) for lo, hi in zip(ends[:-1], ends[1:], strict=True)]
    def hull(values):
        return min(v.lower() for v in values).union(max(v.upper() for v in values))
    t, c, _ = [hull(axis) for axis in zip(*pieces, strict=True)]
    assert arb(".01022") < t < arb(".01043")
    assert arb(".8416") < c < arb(".8765")
    assert c+t*t/d > arb(".8863")
    print("顶耦合整个边际窗的顶一二阶导包围、驻点曲率下端",
          t.str(8), c.str(8), (c+t*t/d).lower().str(14))

    # 五个截止位置固定，只保留权和为一。半个 l1 距离的球是原子间质量转移的凸包。
    y = arb(".8616")
    top_columns = [3*prefactor*y**4*heat(1, y**2*rho/x) for x in sites]
    heavy_columns = [[-prefactor*heat(1, (z+rho*box)/x) for x in sites] for z in masses]
    directions = [(top_columns[i]-top_columns[j], [row[i]-row[j] for row in heavy_columns])
                  for i in range(5) for j in range(5)]
    def least_curvature(moved):
        lows, highs = [], []
        for top_change, heavy_change in directions:
            top = curvature+moved*top_change
            coeff = [w+moved*dw for w,dw in zip(b,heavy_change,strict=True)]
            lows.append(top.lower()+sum(min(arb(0),v.lower()) for v in coeff))
            # 对上端为负的谱行取单位门户，其余取零，即得一个可行上界。
            highs.append(top.upper()+sum(min(arb(0),v.upper()) for v in coeff))
        return min(lows).lower().union(min(highs).upper())
    target = arb(".258")
    left, right = arb(0), arb(2)
    assert least_curvature(left)>target>least_curvature(right)
    for _ in range(24):
        mid = (left+right)/2
        test = least_curvature(mid)
        if test>target:
            left=mid
        elif test<target:
            right=mid
        else:
            break
    assert right-left < arb(".00002")
    assert arb(".77112150") < left < right < arb(".77112199")
    print("仅匹配顶与重标量曲率时归一五权的最小半 l1 改动",
          left.str(13), right.str(13))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
