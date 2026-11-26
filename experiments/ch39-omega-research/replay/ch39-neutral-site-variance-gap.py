"""十九点软辐边球壳：由单点中性方差给两步质量指定共同分界。

依赖同版本 replay/；可显式把其目录作为第一个命令行参数。
"""
from collections import Counter
from fractions import Fraction
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = (Path(sys.argv[1]) if len(sys.argv) > 1 else
        Path(__file__).resolve().parent)
joint = run_path(str(root / "ch39-first-three-mode-observables.py"))
inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))


def neutral_site_variance():
    soft = shell["soft"]
    pole, tuned = soft["pole"], soft["star"]["orbit"]
    degree, softened = Counter(), Counter()
    for i, j, _, _ in shell["geo"]["edges"]:
        degree.update((i, j))
        if pole in (i, j) and (i in tuned or j in tuned):
            softened.update((i, j))
    types = Counter((degree[i], softened[i]) for i in range(19))
    assert types == {(10, 1): 3, (10, 0): 3, (14, 0): 3,
                     (16, 0): 6, (18, 0): 3, (18, 3): 1}
    # 各点方差 = degree-k+k(1-delta)^2；最低型 (10,1) 在 0<delta<1 取到。
    max_delta = arb(23)/25 * arb(2)**-20 / 14**2
    assert soft["point"]["study"]["offset"] > 0 and acb.zeta_zero(1).imag > 14
    assert 9+(1-max_delta)**2 > arb("9.99999999")
    # D 与单位辐边的两步算子在算子范数下的差：L≤36，||ΔL||≤4δ，Δc=δ/3。
    perturbation = ((2*36+15)*4+12) * max_delta
    assert perturbation < arb(".00001")
    print("中性站点（度，软边数）及重数", sorted(types.items()))


def corner_count(q, level):
    fixed_positive, coeff = inertia["schur_coefficients"](q, level)
    radius = q*arb.pi()/16

    @cache
    def phase(n, index):
        # θ 遍历完整闭小格，以球算术同时包围正弦和余弦。
        x = (arb(2*index+1, 1)/n-1)*radius
        return acb(0, x).exp()

    schur_positive = 0
    for constant, a, b, c in coeff:
        center = constant+2*(a+b+c).real
        assert not center.contains(0)
        sign = 1 if center > 0 else -1
        schur_positive += int(sign > 0)
        pending = [(1, 0, 0)]
        while pending:
            n, i, j = pending.pop()
            u, v = phase(n, i), phase(n, j)
            value = constant+2*(a*u+b*v+c*u.conjugate()*v).real
            if sign*value > 0:
                continue
            assert n < 4096, (q, level, sign, i, j, n)
            pending.extend((2*n, 2*i+x, 2*j+y) for x in (0, 1) for y in (0, 1))
    return 19-fixed_positive-schur_positive


def mass_preimage_levels(t):
    lo = (15-(225-4*t).sqrt())/2
    hi = (15+(225+4*t).sqrt())/2
    return lo, 15-lo, hi


def centered_box():
    # 单位辐边时，轻谱由 [0,x₁) 和 (x₂,x₃) 构成。真实质量与其误差 <10⁻⁵。
    for q in (1, 2, 3, 4, 6):
        expected = (0, 7, 10) if q in (1, 3) else (1, 7, 9)
        for t in (arb("8.99998"), arb("10.19002")):
            lo, mid, hi = mass_preimage_levels(t)
            # 荷六最低线可穿过 lo，用 2 包住低段的最大根数，还由 mid、hi 认证恰两根。
            if q == 6:
                assert lo < 2
                lo = arb(2)
            counts = tuple(corner_count(q, x) for x in (lo, mid, hi))
            assert counts == expected, (q, t, counts)
        print("荷、谱分段的三处计数", q, expected)


def actual_witness():
    number = 148341
    phys = joint["joint"]
    model, initial = joint["model_inputs"]()
    ctx.prec = 192
    levels = model["shell"]["shell_spectra"]((0, *model["charges"]), zero_index=number)
    center = sum(levels[0])/18
    ordered = phys["chiral"]["ordered"]
    bands = {q: ordered([v*abs(v-center) for v in values]) for q, values in levels.items()}
    g = acb.zeta_zero(number).imag
    delta = arb(23)/25*arb(2)**-20/(g*g+shell["soft"]["point"]["study"]["offset"])
    cut = 9+(1-delta)**2
    left = max(bands[q][2] for q in model["charges"])
    right = min(bands[q][3] for q in model["charges"])
    assert left < cut < right
    assert bands[6][1] < 9 < bands[6][2] and 9 < left < cut < right
    assert all(bands[q][2] < 9 < bands[q][3] for q in (1, 2, 3, 4))
    assert arb("9.47490790") < left < arb("9.47490792")
    assert arb("10.59169905") < right < arb("10.59169907")
    angles = shell["geo"]["phases"](g)
    difference = abs((angles[724]/angles[804]).arg())
    flavor = phys["flavor"]["singlet_at"](angles)[1]
    assert arb(".23147") < difference < arb(".23148")
    assert arb(".9721712073") < flavor < arb(".9721712074")
    print("两钟角距与六面味坐标", difference.str(14), flavor.str(14))
    gauge = phys["gauge"]
    _, _, alpha = gauge["crossing_enclosure"](model, initial, bands, False, 724,
                                               "29.04", ".1184", recenter=True)
    vac = phys["vacuum"]
    planck = arb("1.22089e19")
    heavy = 32*arb("2.435e18")/724**2*bands[2][15]
    mu2 = (arb("246.22")**2/(heavy*planck))**2
    fields = vac["paired_spectrum"](vac["minimal"]["fields"](bands), min_gap="1e-18")
    rows = fields+vac["low_sector"](4, mu2)
    assert vac["newton"](rows, arb(".03")) < 1 < vac["newton"](rows, arb(".033"))
    coeff = (9*arb(3).log()-12*arb(2).log())/2
    curve = mu2**2*coeff/(arb.pi()*arb("1e-122"))*(1-arb(0).union(mu2/(arb(".03")*coeff)))
    print("零点编号、共同间隙、强耦合与四层曲率比", number, left.str(14), right.str(14),
          alpha.lower().str(11), alpha.upper().str(11), curve.str(14))
    assert arb(".11829042") < alpha < arb(".11829659")
    assert arb(".99835038010") < curve < arb(".99835038011")


def independent_edge_errors():
    for weight_ppm, phase_microradians in ((50, 0), (0, 20), (25, 10)):
        epsilon = Fraction(weight_ppm, 10**6)
        rho = Fraction(phase_microradians, 10**6)
        laplacian_change = 36*epsilon + 108*(1+epsilon)*rho
        mass_change = (87+36*epsilon)*laplacian_change + 540*epsilon*(1+epsilon)
        assert 9 + mass_change < Fraction("9.99999999")*(1-epsilon)**2
        assert 10*(1+epsilon)**2 < Fraction("10.19") - mass_change


def cuts_with_multiple_sites():
    # 取较小一侧，每个顶点在侧内至多有 size-1 个邻点。
    assert min(size*(11-size) for size in range(2, 10)) == 18
    delta_cap = Fraction(23, 25*2**20*14**2)
    assert 18*(1-delta_cap)**2 > 10


if __name__ == "__main__":
    ctx.prec = 192
    cuts_with_multiple_sites()
    independent_edge_errors()
    neutral_site_variance()
    centered_box()
    actual_witness()
