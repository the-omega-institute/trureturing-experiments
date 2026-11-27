"""一阶共同标度驻定的匹配截止测度：总变差下确界的两侧区间证书。

模型与 ch39-cutoff-scale-stationarity.py 相同（三枚候选的可见谱，零裸项，无隐单态，c = 1e-122）。
在 (0,1] 上的有限有符号测度 μ 满足 K_μ=1、W_μ=c、DK_μ=DW_μ=0 当且仅当
    ∫ v dμ = β，v = (1, L, Φ₀, Ψ₀, P)，β = (1, −1, 0, 0, 4c)
（ch39-cutoff-scale-stationarity 说明的定理二）。对任意 y ∈ R⁵，令 F = y·v，则 y·β = ∫ F dμ ≤ sup|F| · ‖μ‖，故
    ‖μ‖ ≥ y·β / sup_{(0,1]} |F|。
y 取自网格线性规划的对偶解（只作为证书，不要求最优）；sup|F| 用 Taylor 区间与二分认证。
上界由显式四点测度给出（用ch39-cutoff-scale-stationarity.py 的 construct：固定三个截止，由行列式变号夹出第四个）。
"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
previous = run_path(str(ROOT / "ch39-cutoff-scale-stationarity.py"))
C_TARGET = arb("1e-122")
DELTA = arb("1e-4")
T_DEEP = arb(-40)

# 网格线性规划的对偶向量（十五位有效数字）与原始解的五个截止 log x。
DUAL = {
    11661: ("-0.999702469910179", "-316340.959488831", "150909.322940951",
            "-0.0134433325308522", "332586598430.911"),
    11934: ("-0.999691544142093", "-328863.110791943", "150512.00475585",
            "-0.0112990586169358", "322715388876.735"),
    12207: ("-0.999702652219035", "-315927.970768828", "150933.506109444",
            "-0.0135160980674422", "333258311991.029"),
}
# 四点上界：三个截止固定（log x = −30 只承担归一），第四个由行列式变号夹出。
FOUR = {
    11661: (("-30", "-1.098", "0"), "-12.45", "-12.35"),
    11934: (("-30", "-1.2425", "0"), "-12.52", "-12.42"),
    12207: (("-30", "-1.0935", "0"), "-12.44", "-12.34"),
}
# 断言区间：对偶下界、四点上界、四点的末位 log x。
EXPECT = {
    11661: (("316308.3", "316308.4"), ("316386.9", "316387.1"), ("-12.394769", "-12.394768")),
    11934: (("328829.2", "328829.3"), ("328912.8", "328912.9"), ("-12.470610", "-12.470609")),
    12207: (("315895.3", "315895.4"), ("315973.9", "315974.0"), ("-12.391939", "-12.391938")),
}
# ch39-cutoff-scale-stationarity 证书的四点解总变差（下界端）。
# 负部质量 b⁻ = (b−1)/2 的下界，由本脚本的对偶下界 b₋ 推出（向下取整）。
NEGATIVE = {11661: "158153", 11934: "164414", 12207: "157947"}
# ch39-cutoff-scale-stationarity.py 四点一阶驻定解的总变差：取其 EXPECT_FIRST 断言区间的下端。
PREVIOUS_FOUR = {n: v[0][0] for n, v in previous["EXPECT_FIRST"].items()}


def decreasing(lo, hi, fn):
    return fn(hi).union(fn(lo))


def heat(rows, tlo, thi):
    """t ∈ [tlo,thi] 上的 Φ_m、Ψ_m（m ≤ 2）、L、P、Q = (1/2π)Σ g z² E₁(u) 的包络。"""
    pi = arb.pi()
    xlo, xhi = tlo.exp(), thi.exp()
    x = xlo.union(xhi)
    ph, ps = [arb(0)] * 3, [arb(0)] * 3
    L = P = Q = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            ph[0] += a * x
            ps[0] += g * x * x
            continue
        ulo, uhi = arb((z / xhi).lower()), arb((z / xlo).upper())
        e = decreasing(ulo, uhi, lambda v: (-v).exp())
        e1 = decreasing(ulo, uhi, lambda v: v.expint(1))
        e2 = decreasing(ulo, uhi, lambda v: v.expint(2))
        for m in range(3):
            ph[m] += a * z**m * e * x**(1 - m) if m < 2 else a * z * z * e / x
            ps[m] += g * z**m * e * x**(2 - m)
        L += a * z * e1
        P += g * z * e2
        Q += g * z * z * e1
    return ph, ps, L / (2 * pi), x * P / (2 * pi), Q / (2 * pi)


def jets(y, rows, tlo, thi):
    """F、F′、F″（′ = d/d log x），用 DL=Φ₁/2π，DΦ₀=Φ₀+Φ₁，DΨ₀=2Ψ₀+Ψ₁，DP=P+Q，DQ=Ψ₂/2π。"""
    pi = arb.pi()
    ph, ps, L, P, Q = heat(rows, tlo, thi)
    F = y[0] + y[1] * L + y[2] * ph[0] + y[3] * ps[0] + y[4] * P
    F1 = y[1] * ph[1] / (2 * pi) + y[2] * (ph[0] + ph[1]) + y[3] * (2 * ps[0] + ps[1]) + y[4] * (P + Q)
    F2 = (y[1] * ph[2] / (2 * pi) + y[2] * (ph[0] + ph[1] + ph[2])
          + y[3] * (4 * ps[0] + 3 * ps[1] + ps[2]) + y[4] * (P + Q + ps[2] / (2 * pi)))
    return F, F1, F2


def sup_certificate(y, rows):
    """证明 [−40,0] 上 |F| ≤ 1+δ：|F(t)| ≤ |F(c)| + |F′(c)|r + ½ sup_B|F″| r²。"""
    count = [0]

    def check(lo, hi, depth):
        count[0] += 1
        c, r = (lo + hi) / 2, (hi - lo) / 2
        F, F1, _ = jets(y, rows, c, c)
        _, _, F2 = jets(y, rows, lo, hi)
        bound = abs(F) + abs(F1) * r + abs(F2) * r * r / 2
        if bound < 1 + DELTA:
            return
        assert depth < 16, ("未能认证", lo.str(10), hi.str(10), bound.str(5))
        check(lo, c, depth + 1)
        check(c, hi, depth + 1)

    t, step = T_DEEP, arb(1) / 16
    while t < 0:
        check(t, t + step, 0)
        t = t + step
    return count[0]


def deep_tail(y, rows):
    """x ≤ e^{−40}：|F − y₀| ≤ |y₁||L| + |y₂||Φ₀| + |y₃||Ψ₀| + |y₄||P|，各项按 u ≥ z_min e^{40} 放大。"""
    xd = T_DEEP.exp()
    zmin = min((z.lower() for _, _, z in rows if not z.is_zero()), key=float)
    U = arb(zmin) / xd
    assert U > 1000
    eU = (-U).exp()
    sa = sum(abs(a) for _, a, z in rows if not z.is_zero())
    sg = sum(abs(g) for g, _, z in rows if not z.is_zero())
    pi = arb.pi()
    # u E₁(u) ≤ e^{−u}，u E₂(u) ≤ e^{−u}：|L| ≤ x·Σ|α|e^{−U}/2π，|P| ≤ x²Σ|g|e^{−U}/2π。
    a0 = sum(abs(a) for _, a, z in rows if z.is_zero())
    g0 = sum(abs(g) for g, _, z in rows if z.is_zero())
    bound = (abs(y[0]) + abs(y[1]) * xd * sa * eU / (2 * pi) + abs(y[2]) * xd * (a0 + sa * eU)
             + abs(y[3]) * xd * xd * (g0 + sg * eU) + abs(y[4]) * xd * xd * sg * eU / (2 * pi))
    assert bound < 1
    return bound


def inside(value, bounds):
    lo, hi = (arb(v) for v in bounds)
    return lo < value < hi


def verify():
    data = previous["spectra"]()
    for n, rows in data.items():
        ctx.prec = 128
        y = [arb(v) for v in DUAL[n]]
        tail = deep_tail(y, rows)
        boxes = sup_certificate(y, rows)
        value = y[0] - y[1] + 4 * C_TARGET * y[4]
        lower = value / (1 + DELTA)
        ctx.prec = 512
        fixed, lo, hi = FOUR[n]
        xs, w4, data4, t4 = previous["construct"](rows, [arb(t) for t in fixed], lo, hi, 1)
        b4 = sum(abs(v) for v in w4)
        # 回代核对 Σw=1、K=1、W=c、ΣwΦ₀ ∋ 0、ΣwΨ₀ ∋ 0（后两项即 DK=DW=0）。
        K = sum(wj * previous["heat"](rows, x, 1)[0] for wj, x in zip(w4, xs))
        W = sum(wj * previous["heat"](rows, x, 1)[1] for wj, x in zip(w4, xs))
        assert (sum(w4) - 1).contains(0) and (K - 1).contains(0) and (W - C_TARGET).contains(0)
        assert sum(wj * previous["heat"](rows, x, 1)[2][0] for wj, x in zip(w4, xs)).contains(0)
        assert sum(wj * previous["heat"](rows, x, 1)[3][0] for wj, x in zip(w4, xs)).contains(0)
        e_lo, e4, et = EXPECT[n]
        assert inside(lower, e_lo) and inside(b4, e4) and inside(t4, et)
        assert lower < b4 and (b4 - lower) / lower < arb("2.6e-4")
        # 与 ch39-cutoff-scale-stationarity 比较：总变差下界 8194；其四点解的总变差 ≥ 4.69e5，是下确界的 1.48 倍以上。
        assert lower > 38 * arb(8194)
        ratio = arb(PREVIOUS_FOUR[n]) / b4
        assert ratio > arb("1.48")
        # 负部质量 b⁻ = (b−1)/2 的下界；四点权的量级。
        assert (lower - 1) / 2 > arb(NEGATIVE[n])
        assert arb("-1.65e5") < w4[0] < arb("-1.57e5") and arb("1.57e5") < w4[3] < arb("1.65e5")
        assert arb("-2.05") < w4[1] < arb("-2.0") and arb("0.16") < w4[2] < arb("0.23")
        print(n, "尾段 sup|F| ≤", tail.str(8), "；Taylor 区间", boxes, "个，sup|F| ≤ 1+1e-4")
        print("   总变差下确界 ∈ [", lower.lower().str(10), ",", b4.upper().str(10), "]",
              "相对宽度", ((b4 - lower) / lower).upper().str(3))
        print("   四点可行测度：log x =", *fixed, t4.mid().str(10), "权重", *(v.str(8) for v in w4))
        print("   ch39-cutoff-scale-stationarity 四点解总变差 / 本上界 >", ratio.lower().str(4))

if __name__ == "__main__":
    verify()
