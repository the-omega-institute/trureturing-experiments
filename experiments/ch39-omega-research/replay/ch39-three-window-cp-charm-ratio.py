"""三窗零点上平均 |Vub|² 精确时，留出 CP 亏额与平均 |Vcb|² 偏移之比的双侧界（上端锐，下端夹到相对 2.6×10⁻⁵）。

记 U=|Vub|²、C=|Vcb|²、J 为各背景按六面坐标 c 读出的味读数，c₀ 为参考相位。
对任意支撑在三窗零点上的味概率 ν，若 E_ν U=U₀ 且 ν 不全落在 c₀，则

    0.0191919011 < D / (10¹²(E_ν C − C₀)) < 0.0638699834 ,    D = 1 − E_ν J / J₀ ,

上端可逼近到 10⁻¹⁰ 以内，下端被夹到相对 2.6×10⁻⁵ 以内。两端由三窗相位区间
两端的单点比 R(c)=H̃(c)/B(c) 给出；这里 H̃=D 的单点核减去切线项，B=10¹²(C−C₀−s(U−U₀))。
本脚本认证：
  1. 外界：在三窗必要相位区间 [0.97561882, a] 上 H̃−κ₋B≥0、κ₊B−H̃≥0，只在 c₀ 取等；
  2. 内界：低端开窗内 R<0.0191923898（与外界相差不足 4.9×10⁻⁷），高端 R(a)>0.0638699833；
  3. 相位方差计：B/(c−c₀)² 与 H̃/(c−c₀)² 在整个区间上的上下界；
  4. 与已有八点、十一点、全谱固定偏移结果的相容性。
"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
two_window = run_path(str(root / "ch39-three-window-two-probability-cp.py"))
extreme = two_window["previous"]
outer = tuple(map(arb, two_window["outer_interval"]))

# 外界常数（向外取整）与内界目标。
KAPPA_LOW_OUT = arb("0.0191919011")
KAPPA_HIGH_OUT = arb("0.0638699834")
KAPPA_LOW_IN = arb("0.0191923898")   # 低端可达值的上界
KAPPA_HIGH_IN = arb("0.0638699833")  # 高端可达值的下界
BETA_LOW, BETA_HIGH = arb("30318.0619"), arb("30318.2522")
ETA_LOW, ETA_HIGH = arb("581.861"), arb("1936.427")


def reference_jets():
    jets = two_window["observables"]()
    cp = two_window["cp"]
    c0 = cp["flavor"]["singlet"](11661)[1]
    up, excess, loss = jets(c0)
    # 参考点上三项读数的值为零，B 的一阶导数为零（s 的定义）。
    assert up[0].contains(0) and excess[0].contains(0) and excess[1].contains(0)
    assert loss[0].contains(0) and up[1] < 0
    tangent = loss[1]/up[1]

    def reduced(c):
        # 返回 (H̃, B) 的二阶射流：H̃=loss−tangent·up，二者在 c₀ 均二阶消失。
        u, b, l = jets(c)
        return [l[i]-tangent*u[i] for i in range(3)], [b[i] for i in range(3)]

    return c0, reduced, tangent


def certify_nonnegative(c0, combination, label):
    """证明 combination(c)≥0 于外区间，且只在 c₀ 为零。combination 返回二阶射流。"""
    pending = [outer]
    cells = 0
    while pending:
        left, right = pending.pop()
        cells += 1
        hull = left.union(right).union(c0)
        # (i) 以 c₀ 为展开点：f(c)=f₂(ξ)(c−c₀)²，ξ 在 [c,c₀] 内。
        if combination(hull)[2] > 0:
            continue
        # (ii) 不含 c₀ 的格以格心作二阶 Taylor 包围。
        if not left.union(right).overlaps(c0):
            mid = (left+right)/2
            width = right-mid
            point = combination(mid)
            region = combination(left.union(right))
            residual = point[0] + point[1]*arb(0, width) + region[2]*arb(0).union(width*width)
            if residual > 0:
                continue
        assert right-left > arb("1e-13"), label
        middle = (left+right)/2
        pending.extend(((left, middle), (middle, right)))
    return cells


def box_value(function, c):
    mid = arb(c.mid())
    at_mid = function(mid)
    over = function(c)
    return [a[0] + b[1]*(c-mid) for a, b in zip(at_mid, over, strict=True)]


def verify():
    # 0. 三窗相位必要区间：重跑全格角距排除，再核对两端。
    ctx.prec = 160
    extreme["angular_bound"]()
    c0, reduced, tangent = reference_jets()
    ctx.prec, ctx.cap = 384, 3
    flavor = two_window["cp"]["flavor"]
    from flint import acb
    axis, other = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    assert other.contains(0) and abs(other) < arb("1e-50")
    assert axis*arb(extreme["aperture"]).cos()-abs(other) > outer[0]
    assert (axis*axis+other*other).sqrt() < outer[1]
    slope = two_window["observables"]()(c0)[0][1]
    assert arb("-91.27630") < slope < arb("-91.27629")
    assert arb("-0.43840494667") < tangent < arb("-0.43840494666")
    print("参考六面坐标 c₀、切线系数、A′(c₀)、振幅 a", c0.str(20), tangent.str(15), slope.str(10), axis.str(20))

    # 1. 比值的外界。
    lower = certify_nonnegative(
        c0, lambda c: [h-KAPPA_LOW_OUT*b for h, b in zip(*reduced(c), strict=True)], "lower")
    upper = certify_nonnegative(
        c0, lambda c: [KAPPA_HIGH_OUT*b-h for h, b in zip(*reduced(c), strict=True)], "upper")
    assert (lower, upper) == (69, 55)
    print("比值外界 κ₋≥", KAPPA_LOW_OUT.str(12), "κ₊≤", KAPPA_HIGH_OUT.str(12), "格数", lower, upper)

    # 3. 相位方差计：B 与 H̃ 除以 (c−c₀)²。
    def square(c):
        return [(c-c0)**2, 2*(c-c0), arb(1)]

    checks = [
        (lambda c: [b-BETA_LOW*q for b, q in zip(reduced(c)[1], square(c), strict=True)], "beta-"),
        (lambda c: [BETA_HIGH*q-b for b, q in zip(reduced(c)[1], square(c), strict=True)], "beta+"),
        (lambda c: [h-ETA_LOW*q for h, q in zip(reduced(c)[0], square(c), strict=True)], "eta-"),
        (lambda c: [ETA_HIGH*q-h for h, q in zip(reduced(c)[0], square(c), strict=True)], "eta+"),
    ]
    for function, label in checks:
        certify_nonnegative(c0, function, label)
    assert (BETA_HIGH-BETA_LOW)/BETA_LOW < arb("6.3e-6")
    assert arb("3.327") < ETA_HIGH/ETA_LOW < arb("3.329")
    print("B/(c−c₀)² 位于", BETA_LOW.str(9), BETA_HIGH.str(9),
          "；H̃/(c−c₀)² 位于", ETA_LOW.str(6), ETA_HIGH.str(7))

    # 2. 两个三窗开相位窗的内界。
    coords = extreme["open_windows"]()
    ctx.prec, ctx.cap = 384, 3
    ratios = []
    for c in coords:
        h_box, b_box = box_value(reduced, c)
        assert b_box > 0
        ratios.append(h_box/b_box)
    low_inner, high_inner = ratios
    assert coords[0] < c0 < coords[1]
    assert coords[0].rad() < arb("5e-13") and coords[1].rad() < arb("5e-13")
    assert arb("0.97561989655") < coords[0] < arb("0.97561989658")
    assert abs(coords[1]-axis) < arb("5e-13")
    assert KAPPA_LOW_OUT < low_inner and low_inner < KAPPA_LOW_IN
    assert KAPPA_LOW_IN - KAPPA_LOW_OUT < arb("4.9e-7")
    assert (KAPPA_LOW_IN - KAPPA_LOW_OUT)/KAPPA_LOW_OUT < arb("2.6e-5")
    # 高端窗中心两钟同相，c=a；窗内任意小的同心子窗仍含实际零点，故 R(a) 是可达比值的极限。
    h_top, b_top = [r[0] for r in reduced(axis)]
    top = h_top/b_top
    assert KAPPA_HIGH_IN < top < KAPPA_HIGH_OUT
    assert high_inner < KAPPA_HIGH_OUT
    print("低端窗单点比、高端窗单点比、R(a)", low_inner.str(12), high_inner.str(12), top.str(13))

    # 两点混合逼近：低端窗配 c₀ 右侧、 高端窗配 c₀ 左侧的近参考点。
    jets = two_window["observables"]()
    ctx.prec, ctx.cap = 384, 3
    # 这里的近点取精确坐标 c₀±2^-e，示范极限过程；实际零点的存在性由等分布定理给出。
    for c_star, side, target in ((coords[0], 1, low_inner), (coords[1], -1, high_inner)):
        h_star, b_star = box_value(reduced, c_star)
        u_far = box_value(lambda x: jets(x)[:1], c_star)[0]
        for exponent in (12, 20, 28):
            c_near = c0 + side*arb(2)**-exponent
            u_near = jets(c_near)[0][0]
            h_near, b_near = [r[0] for r in reduced(c_near)]
            assert u_near*u_far < 0
            weight = -u_near/(u_far-u_near)
            ratio = (weight*h_star+(1-weight)*h_near)/(weight*b_star+(1-weight)*b_near)
            assert abs(ratio-target) < 4*arb(2)**-exponent
            print("  两点混合：远点窗与 c₀±2^-%d 的比值" % exponent, ratio.str(12))

    # 独立核对：在三个坐标上直接对厄米矩阵作谱分解，比较 U、C、J 与射流读数。
    cp = two_window["cp"]
    matrix = cp["flavor"]["reference_matrix"]()
    phase = cp["two"]["invariant"]["phase"]
    ref, _, j_ref = cp["matrix_readings"](matrix, cp["flavor"]["singlet"](11661)[0])
    data = cp["two"]["spectral_data"]()
    for c in (outer[0], arb(".99"), outer[1]):
        readings, _, j_value = cp["matrix_readings"](matrix, phase(c))
        u_direct = (readings[6]-ref[6])*10**5
        loss_direct = 1-j_value/j_ref
        h_direct = loss_direct - tangent*u_direct
        u_jet, _, l_jet = [v[0] for v in jets(c)]
        assert u_direct.overlaps(u_jet) and loss_direct.overlaps(l_jet)
        assert h_direct.overlaps(reduced(c)[0][0])
        assert readings[7].overlaps(cp["two"]["probabilities"](c, data)[1][0])
        print("  直接谱分解核对 c=%s：10⁵ΔU、单点亏额" % c.str(10), u_direct.str(10), loss_direct.str(10))

    # 4. 一致性对照：以下比值由已有说明中的十进制数字直接算出，不是新的证书。
    eight = arb(".06083301838")/arb("1.90254596927")
    assert KAPPA_LOW_OUT < eight < KAPPA_HIGH_OUT
    fixed_low, fixed_high = arb(".04901178")/2, arb(".09914521")/2
    assert KAPPA_LOW_OUT < fixed_low < fixed_high < KAPPA_HIGH_OUT
    # 旧的全相位倍率 1.54×10⁻¹¹ < ΔC/D < 1.51×10⁻⁹ 被收紧为下列锐常数。
    sharp_low, sharp_high = 1/(KAPPA_HIGH_OUT*10**12), 1/(KAPPA_LOW_OUT*10**12)
    assert arb("1.54e-11") < sharp_low < sharp_high < arb("1.51e-9")
    assert arb("1.56568e-11") < sharp_low and sharp_high < arb("5.21054e-11")
    assert arb("98.0") < arb("1.51e-9")/arb("1.54e-11") < arb("98.1")
    assert arb("28.9") < arb("1.51e-9")/arb("5.21054e-11") < arb("29.0")
    # 4.901178% 是全谱下确界的外下界，故只得 D_min/δ 的下界。
    assert arb("2.4505e10") < fixed_low*10**12
    assert 10**12*KAPPA_LOW_IN < arb("2.4505e10")
    # 已有近最小尺度界在两个方向上的换算（均不改进已有界）。
    assert arb("3.81702e-11")*sharp_low > arb("5.97e-22")
    assert arb("3.81702e-11")*sharp_low < arb("5.98e-22") < arb("1.005e-21")
    assert arb("1.92e-11") < arb("1.005e-21")/sharp_high < arb("3.81702e-11")
    print("八点比值、固定偏移两端比值", eight.str(8), fixed_low.str(8), fixed_high.str(8))
    print("ΔC/D 的界", sharp_low.str(6), sharp_high.str(6))


if __name__ == "__main__":
    verify()
