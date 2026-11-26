"""逐质量层的真空配平、三阶差分和规范候选重谱的质量检验。"""
from collections import defaultdict
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


minimal = run_path(str(Path(__file__).with_name("ch39-chiral-minimal-neutral-gravity.py")))
heat = minimal["expint"]


def paired_spectrum(rows, *, min_gap="3e-11"):
    assert arb(min_gap) > 0
    positive = [z for g, _, z in rows if g > 0]
    negative = [z for g, _, z in rows if g < 0]
    assert min(abs(z - w).lower() for z in positive for w in negative) > arb(min_gap)
    assert all(g % 4 == 0 for g, _, _ in rows if g > 0)
    assert sum(-g for g, _, _ in rows if g < 0) == 958
    assert sum(g // 4 for g, _, _ in rows if g > 0) == 162
    # 权重按半单位记。正权的质量保守地推迟，负权的质量保守地提前。
    residual = []
    events = defaultdict(int)
    for g, a, z in rows:
        extra = arb(-g) / 6 if g < 0 else arb(g) / 12
        weight = int((2 * (a + extra)).unique_fmpz())
        residual.append((arb(weight) / 2, z))
        if z.is_zero():
            assert weight == 21
            continue
        bound = z * 10**15
        # 即使证书包围跨越一个整数，也取外侧端点取整。
        cut = bound.upper().ceil() if weight > 0 else bound.lower().floor()
        events[int(cut.unique_fmpz())] += weight
    total = 0
    for edge in sorted(events):
        total += events[edge]
        assert total >= 0
    assert total + 21 == 393
    return residual


def newton(rows, x):
    return x * sum(a * heat(2, z / x) for a, z in rows) / (2 * arb.pi())


def calibrated_cutoff(rows):
    lo, hi = arb(".03"), arb(".033")
    assert newton(rows, lo) < 1 < newton(rows, hi)
    for _ in range(55):
        mid = (lo + hi) / 2
        value = newton(rows, mid)
        if value < 1:
            lo = mid
        elif value > 1:
            hi = mid
        else:
            break
    return lo.union(hi)


def low_sector(degree, mass_squared):
    fermions = (3 * degree + 3) // 4
    scalars_middle = 4 * fermions - 3 * degree
    # 从原配平的零质量层迁出 degree 个实标量，三层质量平方分别为 a、2a、3a。
    scalar_changes = [-degree, 3 * degree, scalars_middle, degree]
    dirac_changes = [0, 0, fermions, 0]
    weights = [s - 4 * f for s, f in zip(scalar_changes, dirac_changes, strict=True)]
    assert weights == [-degree, 3 * degree, -3 * degree, degree]
    assert all(sum(g * j**k for j, g in enumerate(weights)) == 0 for k in range(3))
    return [(arb(s) / 6 + arb(f) / 3, j * mass_squared)
            for j, (s, f) in enumerate(zip(scalar_changes, dirac_changes, strict=True))]


def scan(rows, balanced, degree, coefficient):
    planck, vev, target = arb("1.22089e19"), arb("246.22"), arb("1e-122")
    required = planck * (4 * arb.pi() * target / (degree * coefficient)).sqrt().sqrt()
    trial = (required / planck)**2
    cutoff = calibrated_cutoff(balanced + low_sector(degree, trial))
    assert trial / (cutoff * coefficient) < arb("2e-59")
    print("净权倍数", degree, "目标质量 meV 与牛顿截止",
          (required * 10**12).str(11), cutoff.sqrt().str(12))
    for divisor in (1, 2):
        candidates = [(i, vev**2 / (divisor * planck * z.sqrt()))
                      for i, (_, _, z) in enumerate(rows) if z > 0]
        least = min(v.lower() for _, v in candidates)
        largest = max(v.upper() for _, v in candidates)
        # 核验整批候选的 K=1 根均在 x∈(.03,.033)，随后使用 0≤I(0)−I(t)≤t。
        assert newton(balanced + low_sector(degree, (least / planck)**2), arb(".03")) < 1
        assert newton(balanced + low_sector(degree, (largest / planck)**2), arb(".033")) > 1
        values = []
        for i, mu in candidates:
            loss = (mu / planck)**2 / (arb(".03") * coefficient)
            assert loss < arb("1e-50")
            ratio = (mu / required)**4 * (1 - arb(0).union(loss))
            values.append((abs(ratio - 1), i, ratio))
        winner = min(values, key=lambda entry: entry[0].lower())
        for diff, i, _ in values:
            if i != winner[1]:
                assert winner[0] < diff
        assert winner[0] > arb(".0086")
        print("质量公式分母系数、最接近的行、目标比值",
              divisor, winner[1], winner[2].str(10))


def verify():
    chiral = minimal["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 320
    read = chiral["source"]["spectral_supplier"](model)
    coefficient = (9 * arb(3).log() - 12 * arb(2).log()) / 2
    # 在中等 u 上独立比较四项 E₃ 有限和与替换变量后的正积分。
    u = arb(1) / 7
    direct = (arb(1) / 2 - 3 * heat(3, u) + 3 * heat(3, 2*u) - heat(3, 3*u)) / u**2
    integral = acb.integral(lambda t, analytic: (-t).hypgeom_1f1(1, 2)**3, 0, u)
    assert (direct + integral.real).overlaps(coefficient)
    assert integral.real.is_finite()
    for number in (11661, 11934, 12207):
        source = read(acb.zeta_zero(number).imag)
        bands = {q: chiral["ordered"]([h for h, charge in source if charge == q])
                 for q in (0, *model["charges"])}
        rows = minimal["fields"](bands)
        balanced = paired_spectrum(rows)
        print("正零点与配平牛顿截止", number, calibrated_cutoff(balanced).sqrt().str(15))
        for degree in (1, 4):
            scan(rows, balanced, degree, coefficient)


if __name__ == "__main__":
    verify()
