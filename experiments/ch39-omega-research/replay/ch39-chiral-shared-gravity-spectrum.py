"""共同引力截止下，三枚规范零点至少有两个隐单态质量须向相反方向改变。"""
from itertools import pairwise
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


fixed = run_path(str(Path(__file__).with_name("ch39-chiral-fixed-cutoff-singlets.py")))
heat = fixed["E"]
cutoff_floor = arb(18) / 100


def bound_cutoff_and_slope(rows):
    # 隐谱逐项 E₂ ≤ 2E₃。两项目标要求 12π/x ≤ 6ΣαE₂ − 2ΣgE₃。
    positive_bound = 6 * sum(a for _, a, _ in rows if a > 0)
    positive_bound += sum(-g for g, _, _ in rows if g < 0)
    first_floor = 12 * arb.pi() / positive_bound
    maximum = max(z.upper() for _, _, z in rows) / first_floor
    upper = sum(6 * a * (1 if a > 0 else heat(2, maximum)) for _, a, _ in rows)
    upper += sum(-2 * g * (arb(1) / 2 if g < 0 else heat(3, maximum)) for g, _, _ in rows)
    assert 12 * arb.pi() / upper > cutoff_floor
    # C=ΛG 与 K=G⁻¹/M_P² 的斜率组合消去全部隐单态。
    maximum = max(z.upper() for _, _, z in rows) / cutoff_floor
    weights = [6 * a - g for g, a, _ in rows]
    response = sum(b * ((-maximum).exp() if b > 0 else 1) for b in weights)
    assert response > 206
    slope = cutoff_floor**2 * response / (2 * arb.pi())
    assert slope > arb("1.06")
    assert slope / (1 + 3 * cutoff_floor + 2 * arb("1e-122")) > arb(".68")


def compare(first, second):
    # 按同一荷与谱序配对。下列差值均按第二枚减第一枚。
    maximum = max(z.upper() for _, _, z in first + second)
    signed_mass = newton_mass = log_mass = arb(0)
    signed_variation = newton_variation = response_variation = arb(0)
    for (g, a, z), (h, b, w) in zip(first, second, strict=True):
        assert g == h and a.contains(b) and b.contains(a)
        if z.is_zero() and w.is_zero():
            continue
        assert z > 0 and w > 0
        signed_mass += g * (w - z)
        signed_variation += abs(g * (w - z))
        newton_mass += a * (w - z)
        newton_variation += abs(a * (w - z))
        response_variation += abs((6 * a - g) * (w - z))
        log_mass += a * (w * (w.log() - 1) - z * (z.log() - 1))

    # 下界覆盖 x ≥ 0.18 的整个半轴。第二个主项对 log(x) 的斜率为 −newton_mass。
    vacuum = signed_mass - (1 - heat(2, maximum / cutoff_floor)) * signed_variation
    assert newton_mass < 0
    newton = log_mass + (arb.const_euler() - cutoff_floor.log()) * newton_mass
    newton -= maximum / cutoff_floor * newton_variation
    response = signed_mass - 6 * newton_mass
    response -= (1 - (-maximum / cutoff_floor).exp()) * response_variation
    return vacuum, newton / (2 * arb.pi()), response / (4 * arb.pi())


def mass_ratio(t):
    """已知 E₃(u)=t，包围质量平方比 u。此处两场的 u 均在 [0,4]。"""
    lo, hi = arb(0), arb(4)
    assert heat(3, hi) < t < heat(3, lo)
    for _ in range(60):
        mid = (lo + hi) / 2
        value = heat(3, mid)
        if value > t:
            lo = mid
        elif value < t:
            hi = mid
        else:
            break
    return lo.union(hi)


def second_moment(t):
    return heat(2, mass_ratio(t))


def two_masses(rows):
    q = arb(3) / 7
    s, t = fixed["required"](rows, q)
    # 309 份在所有背景取同一质量 M=Q/32，余下两份分担两个矩。
    s -= 309 * heat(3, arb(1) / 1024)
    t -= 309 * heat(2, arb(1) / 1024)
    assert arb(1) / 2 < s < 1
    low = 2 * second_moment(s / 2)
    high = 1 + second_moment(s - arb(1) / 2)
    assert low < t < high
    # 凸性给单调根。将两个质量一起传入固定质量的截止检验。
    left, right = s / 2, arb(1) / 2
    for _ in range(35):
        mid = (left + right) / 2
        value = second_moment(mid) + second_moment(s - mid)
        if value < t:
            left = mid
        elif value > t:
            right = mid
        else:
            break
    split = left.union(right)
    hidden = [(309, q**2 / 1024), (1, q**2 * mass_ratio(split)),
              (1, q**2 * mass_ratio(s - split))]
    witness_zero(rows, hidden, q**2)
    return s, t, low, high


def witness_zero(rows, hidden, x0):
    # 前面的通用符号界保证首次过零在这个较宽的区间内。
    x = x0.union(x0 * (1 + arb("1e-6")))
    masses = [(g, z) for g, _, z in rows] + hidden
    moment = sum(g * z * heat(2, z / x) for g, z in masses)
    c = arb("1e-122")
    delta = 4 * arb.pi() * c / (x0 * moment - 4 * arb.pi() * c)
    relative = delta / (1 + (1 + delta).sqrt())
    assert relative > 0 and relative < arb("1e-118")
    print("311 单态实例的首次过零相对截止增量", relative.str(9))


def positive_curvature_window(rows):
    # x0 是某次同时匹配 K 与 C 的截止。上一界给 x0 > 0.18。
    maximum = max(z.upper() for _, _, z in rows) / cutoff_floor
    e2, e3 = heat(2, maximum), heat(3, maximum)
    a_min = sum(a * (e2 if a > 0 else 1) for _, a, _ in rows)
    a_max = sum(a * (1 if a > 0 else e2) for _, a, _ in rows)
    g_max = sum(g * (arb(1) / 2 if g > 0 else e3) for g, _, _ in rows)
    c = arb("1e-122")
    # 隐扇区 T>S 故有上界，与单态数量、质量分布无关。
    denom = 6 * a_min - g_max - 4 * arb.pi() * c / cutoff_floor**2
    assert denom > 0
    x_max = 12 * arb.pi() / denom
    assert x_max < arb(".726")
    exp_min = sum(a * ((-maximum).exp() if a > 0 else 1) for _, a, _ in rows)
    # x≥x0 时，新增单态的 exp(-z/x) 至少为其在 x0 处的 E₂。
    newton_slope = (2 * arb.pi() / x_max - a_max + exp_min) / (2 * arb.pi())
    assert newton_slope > arb("1.37")
    mass_moment = sum(g * z * (e2 if g > 0 else 1) for g, _, z in rows)
    assert mass_moment > arb("4.88e-5")
    # 用 δ=r(2+r) 表示 x/x0−1，避免在 1+r 处损失 117 位小数。
    r = arb("7.2e-117")
    delta = r * (2 + r)
    assert c * (1 + delta) < mass_moment * cutoff_floor * delta / (4 * arb.pi())
    print("可行 x 上界、K′ 下界与有符号质量矩下界",
          x_max.upper().str(9), newton_slope.lower().str(9), mass_moment.lower().str(9))


def verify():
    chiral = fixed["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    spectra = []
    for number in (12207, 11661, 11934):
        source = read(acb.zeta_zero(number).imag)
        bands = {q: chiral["ordered"]([h for h, charge in source if charge == q])
                 for q in (0, *model["charges"])}
        rows = fixed["minimal"]["fields"](bands)
        bound_cutoff_and_slope(rows)
        positive_curvature_window(rows)
        print(number, "剩余真空矩、牛顿矩与牛顿矩两端",
              *(v.str(10) for v in two_masses(rows)))
        spectra.append((number, rows))
    for (i, left), (j, right) in pairwise(spectra):
        vacuum, newton, response = compare(left, right)
        assert vacuum > arb("9.0e-8") and newton > arb("1.47e-8") and response > arb("1.08e-8")
        print(i, j, "xΔS、可见牛顿倒数差、消去单态的斜率组合差下界",
              *(v.lower().str(9) for v in (vacuum, newton, response)))


if __name__ == "__main__":
    verify()
