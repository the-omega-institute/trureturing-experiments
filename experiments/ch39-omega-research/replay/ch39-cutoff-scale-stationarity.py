"""有符号归一截止测度下牛顿与曲率读数对共同截止标度的驻定：匹配 K=1 时两点不能驻定真空项，
四点给一阶、六点给二阶共同标度驻定；单个截止独立移动时敏感度并未改善。

模型与第十一章相同：电弱对称原点的十六重标量可见谱（#11661、#11934、#12207），
零裸项、无任何隐单态。截止测度 μ = Σ_j w_j δ_{x_j}，Σ w_j = 1，x_j = Q_j²/M_P² ∈ (0,1]，
K_μ(λ) = Σ w_j K(λx_j)，W_μ(λ) = Σ w_j W(λx_j)，C = W/K²，D = λ d/dλ。
"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
minimal = run_path(str(ROOT / "ch39-chiral-minimal-neutral-gravity.py"))
expint = minimal["expint"]
C_TARGET = arb("1e-122")
BOOK_SINGLE_CUTOFF = arb("7.2e-117")


def spectra():
    chiral = minimal["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    out = {}
    for n in (11661, 11934, 12207):
        rows = read(acb.zeta_zero(n).imag)
        layers = {q: chiral["ordered"]([h for h, z in rows if z == q])
                  for q in (0, *model["charges"])}
        out[n] = minimal["fields"](layers)
    return out


def heat(rows, x, order):
    """K, W, Φ_m = x^{1-m}Σαz^m e^{-z/x}, Ψ_m = x^{2-m}Σgz^m e^{-z/x}（m<order）。

    DK = Φ_0/2π，DW = −Ψ_0/4π，DΦ_m = (1−m)Φ_m + Φ_{m+1}，DΨ_m = (2−m)Ψ_m + Ψ_{m+1}。
    """
    pi = arb.pi()
    K = W = arb(0)
    P = [arb(0)] * order
    Q = [arb(0)] * order
    for g, a, z in rows:
        if z.is_zero():
            K += a
            W += arb(g) / 2
            P[0] += a
            Q[0] += g
            continue
        u = z / x
        K += a * expint(2, u)
        W += g * expint(3, u)
        e = (-u).exp()
        for m in range(order):
            P[m] += a * z**m * e
            Q[m] += g * z**m * e
    phi = [x**(1 - m) * P[m] for m in range(order)]
    psi = [x**(2 - m) * Q[m] for m in range(order)]
    return x * K / (2 * pi), -x * x * W / (4 * pi), phi, psi


def two_cutoff_obstruction(rows):
    """两点 μ 不能同时给 K_μ=1 与 DW_μ=0（与 C 取何值无关）。"""
    pi = arb.pi()
    # (i) 按质量的任意前缀净偏振 ≤ −62：区间质量取最坏次序。
    ends = [arb(0)] + [z.upper() for *_, z in rows] + [z.lower() for *_, z in rows]
    worst = max(sum(g for g, _, z in rows
                    if (g > 0 and z.lower() <= t) or (g < 0 and z.upper() <= t))
                for t in ends)
    assert worst == -62
    # (ii) 0<x≤1 时 L=(1/2π)ΣαzE₁(z/x) 满足 |L| ≤ (1/2π)Σ|α|zE₁(z)。
    l_max = sum(abs(a) * z * expint(1, z) for _, a, z in rows if not z.is_zero()) / (2 * pi)
    assert l_max < arb("1.2205e-4")
    # (iii) K<1：x<x_c 时 K<(x/2π)Σ_{α>0}α<1；x≥x_c 时 ΣαE₂(z/x) 的上界已为负。
    positive = sum(a for _, a, _ in rows if a > 0)
    x_c = 2 * pi / positive * arb("0.99")
    upper = positive + sum(a * expint(2, z / x_c) for _, a, z in rows if a < 0)
    assert upper < arb("-17.15")
    # 按 (i)，Σg e^{-u} ≤ −62 且 Σg(1+u)e^{-u} ≤ −62，于是 DΨ₀/Ψ₀ > 1；
    # 由 1−K = 1+L−DK 与 (ii)(iii)，D log[Ψ₀/(K−1)] > 0：比值严格单调，两点必重合。
    return l_max


def first_order_cost_floor(l_max, bounds):
    # K_μ=1 且 DK_μ=0 ⇔ Σw L(x_j) = −1，故总变差 b ≥ 1/sup|L|（此界不锐）。
    floor = 1 / l_max
    assert floor > arb(bounds)
    return floor


def two_point_without_newton(rows):
    """不要求 K=1 时两点可以同时驻定：Φ₀/Ψ₀ 在 log x≈−13 附近转向。"""
    def ratio(t):
        _, _, phi, psi = heat(rows, t.exp(), 1)
        return phi[0] / psi[0], phi[0], psi[0]
    t2 = arb("-12.5")
    target, phi2, psi2 = ratio(t2)
    lo, hi = arb("-13.5"), arb("-13.25")
    f_lo, f_hi = ratio(lo)[0] - target, ratio(hi)[0] - target
    assert (f_lo < 0 and f_hi > 0)
    for _ in range(80):
        mid = (lo + hi) / 2
        f = ratio(mid)[0] - target
        if f.contains(0):
            break
        if f < 0:
            lo = mid
        else:
            hi = mid
    assert ratio(lo)[0] < target < ratio(hi)[0]
    t1 = lo.union(hi)
    x1, x2 = t1.exp(), t2.exp()
    K1, _, phi1, psi1 = heat(rows, x1, 1)
    K2, _, _, _ = heat(rows, x2, 1)
    w2 = phi1[0] / (phi1[0] - phi2)
    w1 = 1 - w2
    # 在比值相等的根处，Σw=1、ΣwΦ₀=0 推出 ΣwΨ₀=0，即 DK_μ=DW_μ=0、DC_μ=0。
    assert (w1 * phi1[0] + w2 * phi2).contains(0) and (w1 * psi1[0] + w2 * psi2).contains(0)
    newton = w1 * K1 + w2 * K2
    assert arb("1.2845") < w1 < arb("1.2846") and arb("2.98e-6") < newton < arb("3.00e-6")
    return t1, w1, newton


def single_position_flip(rows, xs, w, data):
    """只移动第 j 个截止：∂W_μ/∂log x_j = w_j DW(x_j)≠0，求曲率翻负所需的相对 Q 增量上界。"""
    pi = arb.pi()
    eta = arb("1e-3")
    out = []
    for wj, x, d in zip(w, xs, data):
        slope = -wj * d[3][0] / (4 * pi)
        assert not slope.contains(0)
        delta = C_TARGET / abs(slope)
        band = (-2 * delta).exp().union((2 * delta).exp())
        _, _, _, psi = heat(rows, band * x, 2)
        M2 = abs(wj) * abs(2 * psi[0] + psi[1]) / (4 * pi)
        outer = delta * (1 + eta)
        # W(δ) = c + slope·δ + R，|R| ≤ M2 δ²/2；沿使 W 下降的方向 δ=−sign(slope)·outer 已为负。
        assert C_TARGET - abs(slope) * outer + M2 * (2 * delta)**2 / 2 < 0
        assert C_TARGET - abs(slope) * delta * (1 - eta) - M2 * (2 * delta)**2 / 2 > 0
        r = max(abs((outer / 2).exp() - 1), abs(1 - (-outer / 2).exp()))
        out.append(r)
    return out


def second_order_cost_coefficient(rows):
    # DK=D²K=0 ⇒ Σ_j w_j Σαz e^{-z/x_j}=0 ⇒ |Str αz| ≤ b Σ|α|z²/x_min。
    moment = sum(a * z for _, a, z in rows)
    spread = sum(abs(a) * z * z for _, a, z in rows)
    ratio = abs(moment) / spread
    assert ratio > 7.2e5
    return ratio


def vector(rows, x, order):
    K, W, phi, psi = heat(rows, x, order)
    row = [K - 1, W - C_TARGET]
    for m in range(order):
        row += [phi[m], psi[m]]
    return row


def construct(rows, fixed, lo, hi, order):
    """固定 2order+1 个 t_j=log x_j，在 [lo,hi] 中由行列式变号夹出末位截止。"""
    base = [vector(rows, arb(t).exp(), order) for t in fixed]
    size = len(base) + 1

    def det(t):
        last = vector(rows, t.exp(), order)
        cols = base + [last]
        return arb_mat([[cols[j][i] for j in range(size)] for i in range(size)]).det()

    lo, hi = arb(lo), arb(hi)
    d_lo, d_hi = det(lo), det(hi)
    assert (d_lo > 0 and d_hi < 0) or (d_lo < 0 and d_hi > 0)
    for _ in range(64):
        mid = (lo + hi) / 2
        d_mid = det(mid)
        if (d_mid > 0) == (d_lo > 0) and not d_mid.contains(0):
            lo = mid
        elif not d_mid.contains(0):
            hi = mid
        else:
            break
    assert (det(lo) > 0) != (det(hi) > 0)
    t_last = lo.union(hi)
    xs = [arb(t).exp() for t in fixed] + [t_last.exp()]
    data = [heat(rows, x, order + 2) for x in xs]
    vecs = base + [vector(rows, xs[-1], order)]
    # 首行归一，余下 size−1 行取前 size−1 个分量；末分量由奇异性自动满足。
    A = arb_mat([[arb(1)] * size] + [[vecs[j][i] for j in range(size)] for i in range(size - 1)])
    rhs = arb_mat([[arb(1)]] + [[arb(0)] for _ in range(size - 1)])
    w = A.solve(rhs).entries()
    last_component = sum(w[j] * vecs[j][size - 1] for j in range(size))
    assert last_component.contains(0)
    return xs, w, data, t_last


def taylor_bound(rows, xs, w, k, eps0):
    """sup_{|log λ|≤eps0} |D^k W_μ|，用 D^kW = −(1/4π) D^{k−1}Ψ₀ 的展开系数。"""
    coeff = {3: (4, 3, 1, 0), 4: (8, 7, 3, 1)}[k]
    lam = (-eps0).exp().union(eps0.exp())
    total = arb(0)
    for wj, x in zip(w, xs):
        _, _, _, psi = heat(rows, lam * x, 4)
        total += abs(wj) * abs(sum(c * p for c, p in zip(coeff, psi)))
    return total / (4 * arb.pi())


def newton_positive(rows, xs, w, eps0):
    """K(0)=1、DK(0)=0，故 K(ε) ≥ 1 − sup|D²K| ε²/2；窗内曲率与 W 同号。"""
    lam = (-eps0).exp().union(eps0.exp())
    total = arb(0)
    for wj, x in zip(w, xs):
        _, _, phi, _ = heat(rows, lam * x, 2)
        total += abs(wj) * abs(phi[0] + phi[1])
    assert 1 - total / (2 * arb.pi()) * eps0**2 / 2 > 0


def first_order_window(rows, xs, w, data):
    pi = arb.pi()
    # 恒等式核对：K−DK = −L ⇒ Σ w L = −1；W − DW/2 = P/4 ⇒ Σ w P = 4c。
    L = sum(wj * sum(a * z * expint(1, z / x) for _, a, z in rows if not z.is_zero())
            for wj, x in zip(w, xs)) / (2 * pi)
    assert (L + 1).contains(0)
    K2 = sum(wj * d[2][1] for wj, d in zip(w, data)) / (2 * pi)
    W2 = -sum(wj * d[3][1] for wj, d in zip(w, data)) / (4 * pi)
    assert W2 < 0
    eps = (2 * C_TARGET / -W2).sqrt()
    eps0 = 2 * eps
    M3 = taylor_bound(rows, xs, w, 3, eps0)
    eta = arb("1e-4")
    inner, outer = eps * (1 - eta), eps * (1 + eta)
    # W(ε) = c + W2 ε²/2 + R，|R| ≤ M3|ε|³/6（W(0)=c、DW(0)=0 精确成立）。
    assert C_TARGET + W2 * inner**2 / 2 - M3 * eps0**3 / 6 > 0
    assert C_TARGET + W2 * outer**2 / 2 + M3 * eps0**3 / 6 < 0
    newton_positive(rows, xs, w, eps0)
    b = sum(abs(x) for x in w)
    r_inner, r_outer = (inner / 2).exp() - 1, (outer / 2).exp() - 1
    return b, W2, K2, r_inner, r_outer, 1 - (-outer / 2).exp()


def second_order_window(rows, xs, w, data):
    pi = arb.pi()
    W3 = -sum(wj * d[3][2] for wj, d in zip(w, data)) / (4 * pi)
    K3 = sum(wj * d[2][2] for wj, d in zip(w, data)) / (2 * pi)
    assert W3 > 0
    eps = (6 * C_TARGET / W3) ** (arb(1) / 3)
    eps0 = 2 * eps
    M4 = taylor_bound(rows, xs, w, 4, eps0)
    eta = arb("3e-3")
    inner, outer = eps * (1 - eta), eps * (1 + eta)
    # W(ε) = c + W3 ε³/6 + R，|R| ≤ M4 ε⁴/24；ε<0 一侧先翻负，ε>0 一侧在 eps0 内恒正。
    assert C_TARGET - W3 * inner**3 / 6 - M4 * eps0**4 / 24 > 0
    assert C_TARGET - W3 * outer**3 / 6 + M4 * eps0**4 / 24 < 0
    assert C_TARGET - M4 * eps0**4 / 24 > 0
    newton_positive(rows, xs, w, eps0)
    b = sum(abs(x) for x in w)
    r_up = (eps0 / 2).exp() - 1
    assert r_up > arb("1.4e-37")
    return b, W3, K3, 1 - (-inner / 2).exp(), 1 - (-outer / 2).exp(), r_up


FIRST = {
    11661: ((-16, -13, -8), -7, -6),
    11934: ((-18, -12, -11), -9, -8),
    12207: ((-16, -13, -8), -6, -5),
}
SECOND = {
    11661: ((-28, -26, -24, -20, -18), -14, -12),
}
# 总变差、D²W（或 D³W）与翻号相对截止增量的断言区间。
EXPECT_FIRST = {
    11661: (("469698", "469701"), ("-4.8125e-7", "-4.8123e-7"), ("1.0191e-58", "1.0195e-58")),
    11934: (("507233", "507236"), ("-5.1077e-7", "-5.1075e-7"), ("9.892e-59", "9.896e-59")),
    12207: (("469042", "469045"), ("-4.8031e-7", "-4.8029e-7"), ("1.0201e-58", "1.0205e-58")),
}
EXPECT_SECOND = {
    11661: (("1.2747e13", "1.2748e13"), ("2.044e-11", "2.046e-11"), ("7.1355e-38", "7.1796e-38"),
            ("7.4388", "7.4390"), ("-12.23390", "-12.23389")),
}
# D²K 与末位 log x 的断言区间；一阶代价下界。
EXPECT_EXTRA = {
    11661: (("-0.64612", "-0.64610"), ("-6.08730", "-6.08729"), "8195.17"),
    11934: (("-0.56790", "-0.56789"), ("-8.50985", "-8.50984"), "8210.41"),
    12207: (("-0.64533", "-0.64532"), ("-5.30658", "-5.30657"), "8194.12"),
}
# 重谱上端（GeV）与最高截止倍数的断言区间。
EXPECT_RANGE = {
    11661: ("1.266e16", "1.268e16", "45.9", "46.0"),
    11934: ("1.242e16", "1.243e16", "13.9", "14.0"),
    12207: ("1.268e16", "1.270e16", "67.7", "67.8"),
}
# 只移动一个截止时翻负所需相对 Q 增量的最小值上界。
EXPECT_SINGLE = {11661: "1.8e-120", 11934: "8.6e-120", 12207: "2.0e-120"}
EXPECT_SINGLE_SIX = ("3.39e-116", "3.40e-116", "3.43e-112")


def inside(value, bounds):
    lo, hi = (arb(t) for t in bounds)
    return lo < value < hi


def verify():
    data = spectra()
    ctx.prec = 512
    for n, rows in data.items():
        assert sum(g for g, _, _ in rows) == -310
        assert (sum(a for _, a, _ in rows) + arb(103) / 6).contains(0)
        l_max = two_cutoff_obstruction(rows)
        k2_bounds, t_bounds, floor_bound = EXPECT_EXTRA[n]
        floor = first_order_cost_floor(l_max, floor_bound)
        # 仅数值：log x∈[−40,0] 步长 0.1 的格点上 |L| 的最大值（格点估计，不是证书）。
        grid = max(abs(sum(a * z * expint(1, z / arb(-40 + i / 10).exp())
                           for _, a, z in rows if not z.is_zero())).mid() for i in range(401))
        grid = arb(grid) / (2 * arb.pi())
        assert arb("5.2e-5") < grid < arb("5.4e-5")
        print("   仅数值：格点上 max|L| 与相应总变差估计", grid.str(4), (1 / grid).str(4))
        ratio = second_order_cost_coefficient(rows)
        print(n, "两点匹配 K=1 时不能驻定真空项；sup|L| 上界、一阶驻定总变差下界、二阶驻定系数 |Strαz|/Σ|α|z²",
              l_max.upper().str(6), floor.lower().str(6), ratio.lower().str(6))
        fixed, lo, hi = FIRST[n]
        xs, w, heat_data, t_last = construct(rows, fixed, lo, hi, 1)
        b, W2, K2, r_in, r_out, r_down = first_order_window(rows, xs, w, heat_data)
        assert b > floor and floor > 8194 and ratio > arb("7.2e5")
        (bb, ww, rr) = EXPECT_FIRST[n]
        assert inside(b, bb) and inside(W2, ww)
        assert arb(rr[0]) < r_in and r_out < arb(rr[1]) and r_down < arb(rr[1])
        assert r_in / BOOK_SINGLE_CUTOFF > arb("1.374e58")
        assert inside(K2, k2_bounds) and inside(t_last, t_bounds)
        singles = single_position_flip(rows, xs, w, heat_data)
        assert min(singles) < arb(EXPECT_SINGLE[n])
        assert min(singles) < BOOK_SINGLE_CUTOFF / 800
        gev = [x.sqrt() * arb("1.22089e19") for x in xs]
        assert arb("1.4e15") < min(gev) and max(gev) < arb("8.7e17")
        heavy = [z for *_, z in rows if not z.is_zero()]
        assert arb("1.5e-11") < min(heavy) and max(heavy) < arb("1.081e-6")
        top = max(heavy).sqrt() * arb("1.22089e19")
        print("   重谱下端、上端 GeV", (min(heavy).sqrt() * arb("1.22089e19")).str(4), top.str(5), "最高截止倍数", (max(gev) / top).str(4))
        bottom = min(heavy).sqrt() * arb("1.22089e19")
        top_lo, top_hi, fac_lo, fac_hi = EXPECT_RANGE[n]
        assert arb("4.9e13") < bottom < arb("4.94e13")
        assert arb(top_lo) < top < arb(top_hi) and arb(fac_lo) < max(gev) / top < arb(fac_hi)
        # 逐点贡献：w_jK(x_j) 为 O(1)，w_jW(x_j) 相消到 c。
        kc = [wj * heat(rows, x, 1)[0] for wj, x in zip(w, xs)]
        wc = [wj * heat(rows, x, 1)[1] for wj, x in zip(w, xs)]
        kmax = max(abs(v).upper() for v in kc)
        wmax = max(abs(v).upper() for v in wc)
        assert arb("0.84") < kmax < arb("2.6") and arb("2.9e-4") < wmax < arb("1.5e-3")
        assert arb("6.9e-120") < C_TARGET / wmax < arb("3.4e-119")
        print("   逐点 w_jK(x_j)", *(v.str(3) for v in kc), "逐点 w_jW(x_j)", *(v.str(3) for v in wc),
              "c/max|w_jW|", (C_TARGET / wmax).str(3))
        assert arb("2.3e5") < max(abs(x) for x in w) < arb("2.6e5")
        # 平行工作的负部质量 b⁻ 与总变差的关系：b = 1 + 2b⁻。
        assert (floor - 1) / 2 > arb("4096.5")
        print("   四个截止 Q/GeV", *(q.str(4) for q in gev))
        print(n, "四点一阶驻定：末位 log x", t_last.mid().str(12),
              "总变差", b.str(8), "D²W", W2.str(8), "D²K", K2.str(8))
        print("   曲率保持为正的相对截止上调 >", r_in.lower().str(8),
              "，上调", r_out.upper().str(8), "或下调", r_down.upper().str(8), "已为负")
        print("   沿共同标度相对单截止通用界 7.2e-117 的放大倍数 >", (r_in / BOOK_SINGLE_CUTOFF).lower().str(6))
        print("   只移动单个截止时翻负所需相对 Q 增量（上界，逐点）", *(r.upper().str(4) for r in singles))
        if n == 11661:
            t1, w1, newton = two_point_without_newton(rows)
            print("   不要求 K=1 的两点驻定反例：log x₁", t1.mid().str(10), "log x₂ −12.5，w₁",
                  w1.str(8), "K_μ", newton.str(6))
        if n in SECOND:
            fixed, lo, hi = SECOND[n]
            xs, w, heat_data, t_last = construct(rows, fixed, lo, hi, 2)
            b, W3, K3, r_in, r_out, r_up = second_order_window(rows, xs, w, heat_data)
            (bb, ww, rr) = EXPECT_SECOND[n][:3]
            assert inside(b, bb) and inside(W3, ww)
            assert arb(rr[0]) < r_in and r_out < arb(rr[1])
            assert inside(K3, EXPECT_SECOND[n][3]) and inside(t_last, EXPECT_SECOND[n][4])
            singles = single_position_flip(rows, xs, w, heat_data)
            lo6, mid6, hi6 = (arb(v) for v in EXPECT_SINGLE_SIX)
            assert all(r < hi6 for r in singles) and lo6 < min(singles) < mid6
            print(n, "六点二阶驻定：末位 log x", t_last.mid().str(12),
                  "总变差", b.str(8), "D³W", W3.str(8), "D³K", K3.str(8))
            print("   下调", r_in.lower().str(8), "以内曲率为正，下调",
                  r_out.upper().str(8), "已为负；上调", r_up.lower().str(6), "以内恒正")
            print("   只移动单个截止时翻负所需相对 Q 增量（上界，逐点）", *(r.upper().str(4) for r in singles))


if __name__ == "__main__":
    verify()
