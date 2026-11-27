"""第七壳三门引力和弱差共同根的两个规范留出量，并核查中途动能。"""

from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def tighten(readings, weights, lib, p0):
    # 只剔除无法满足 D=0 的区间，保留原三门盒中的每一个精确根。
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    err = sum(
        abs(weights[i]) * (arb("2.398e-21") + arb("7.468e-36") / y**2)
        for i, y in enumerate(readings)
        if i != 15
    )
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    todo = [(readings[15].lower(), readings[15].upper())]
    kept = []
    while todo:
        l, r = todo.pop()
        w = weak(l.union(r) ** 2 / p0)
        if not (
            reference + weights[15] * (w - reference) + arb(0, err.upper())
        ).contains(0):
            continue
        if r - l < arb("1e-25"):
            kept.append(l.union(r))
        else:
            mid = (l + r) / 2
            todo.extend(((l, mid), (mid, r)))
    assert kept
    result = kept[0]
    for u in kept[1:]:
        result = result.union(u)
    assert arb("7.8189e-17") < result < arb("7.8192e-17")
    readings[15] = result
    print("弱差方程收紧后的第十六层振幅", result.str(12), flush=True)


def check_points(root, readings, weights, lib):
    gauge = lib("all-field-signed-crossing")
    num = 2_306_888_100_092
    data = gauge["joint_system"](root, "1.176397e14", number=num, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, num)
    ctx.prec = 192
    neutral = lib("full-neutral-free-cutoff")
    species, _ = neutral["response_rows"](
        data,
        masses,
        signs,
        (arb("246.22") / arb("1.22089e19")) ** 2 / 2,
        arb(".42995401872"),
    )
    tighten(readings, weights, lib, data["p0"])
    sites = [v**2 for v in readings]
    center = [v.mid() for v in sites]
    running = lib("full-neutral-running")
    rule = gauge["angular_rule"](10, 29)
    weak = running["weak_bound"](
        data,
        [data["p0"], *sites],
        [arb(0), *weights],
        lib("neutral-axial-top"),
        lib("neutral-bosonic"),
    )
    variation = sum(
        abs(w) * abs((x / y).log())
        for w, x, y in zip(weights, sites, center, strict=True)
    )
    errs = [
        sum(abs(row[2][j]) for row in species) * variation / (4 * arb.pi())
        for j in range(3)
    ]
    limits = (
        ("3200", -280, -30000),
        ("1e7", -4500, -9_000_000),
        ("1e10", -7900, -28_000_000),
        ("1.176397e14", -12300, -65_000_000),
    )
    for energy, col_upper, det_upper in limits:
        cols = running["direct"](data, center, weights, species, energy, gauge, rule)
        q, m, c = [v + arb(0, e.upper()) for v, e in zip(cols, errs, strict=True)]
        d = arb(0) if energy == "1.176397e14" else arb(0, weak.upper())
        det = m * (q - m) + q * d
        assert c < col_upper and det < det_upper
        if energy == "1.176397e14":
            assert 39900 < q - 8 * c / 3 < 40100 and 6900 < m - c < 7100
        print(
            "能量、QQ、QT3、色、中性行列式、RP、RM",
            energy,
            *[v.str(13) for v in (q, m, c, det, q - 8 * c / 3, m - c)],
            flush=True,
        )


def check_components(root, readings, frame, lib):
    target = lib("fifth-third-pivot-em-colour-path")["solve"](
        root,
        lib("seventh-third-pivot-positive-path")["CENTERS"],
        (17, 0, 15, 1, 4, 6),
        number=2_306_888_100_092,
    )[4]
    floor_amp, _, highest, out, _, _, phases = frame
    cluster = (1, 3, 7, 10, 13, 16)
    # 重相位相同的六行几乎只含一条重方向；对所有酉变化和二行压缩仍有低方向。
    secondary = (highest - floor_amp) * sum(
        abs(out[i, 1] - phases[0] * out[i, 0]) ** 2 for i in cluster
    )
    assert all(
        floor_amp + secondary < target[i].sqrt() < readings[i] for i in cluster[1:]
    )
    print("六层压缩超底的第二谱值上界", secondary.upper().str(11), flush=True)
    bright = (1, 4, 6)
    low_targets = [v.sqrt() for i, v in enumerate(target) if i not in bright]
    leak = sum(v - floor_amp for v in low_targets)
    reservoir = (readings[4] + readings[6] - sum(low_targets)) / 2
    assert all(floor_amp < v < reservoir for v in low_targets)
    largest = max(v.sqrt().upper() for v in target)
    assert highest - leak > largest
    print(
        "逐层设置前的谱上界下限、目标十六层总超底振幅、三层最大谱与目标最大值的间隔",
        reservoir.lower().str(10),
        leak.upper().str(10),
        (highest - leak - largest).lower().str(10),
        flush=True,
    )
    changes = [t.sqrt() - s for t, s in zip(target, readings, strict=True)]
    positive = {i for i, v in enumerate(changes) if v > 0}
    negative = {i for i, v in enumerate(changes) if v < 0}
    assert positive == {1, 4, 6} and len(negative) == 16
    # 若两个以上分量各自守迹，至少一块及其补块都含正差和负差。
    # 每个非平凡二分中，含三个正差中恰好一个的一侧唯一。
    floor = arb(1)
    count = 0
    negative = sorted(negative)
    for size in range(17):
        for group in combinations(negative, size):
            dec = sum((changes[i] for i in group), arb(0))
            for i in positive:
                gap = abs(changes[i] + dec)
                floor = min(floor, gap.lower())
                count += 1
    assert count == 196608 and floor > arb(".01088")
    print("含单个正差的割数、割迹绝对差统一下界", count, floor.str(14), flush=True)


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    source = lib("fifth-seventh-three-gates-weak-gravity")
    data = source["inputs"](root)
    readings = source["three_gates"](root, data)
    check_points(root, readings, data[2], lib)
    check_components(root, readings, data[-1], lib)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
