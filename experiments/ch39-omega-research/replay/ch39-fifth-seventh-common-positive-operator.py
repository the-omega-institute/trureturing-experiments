"""同一个十九维荷六正算子不能同时产生第五壳全程正解与第七壳真空目标。"""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root, *, sixth_measurement=False, sixth_schur=False, readout_radius=0, swap_five_sixteen=False):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    fifth, seventh = 2_306_888_100_090, 2_306_888_100_092
    source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2][1]
    weights = []
    native = lib("native-rh-ray-pvm")
    for number in (fifth, seventh):
        pivots = next(v for v in native["shell_weights"](root, number) if v[0] == 6)[2]
        weights.append([v.real for v in pivots[2]])
    c, e = weights
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    bases = [spectrum(shell, 6, n)[1] for n in (fifth, seventh)]
    overlap = bases[1].conjugate().transpose() * bases[0]
    if sixth_measurement or sixth_schur:
        middle = spectrum(shell, 6, fifth + 1)[1]
        incoming = middle.conjugate().transpose() * bases[0]
        outgoing = bases[1].conjugate().transpose() * middle
    _, heat, fields = lib("three-index-free-atoms")["broken_fields"](
        root, [fifth, seventh]
    )
    ctx.prec = 192
    radius = arb(readout_radius)
    assert radius >= 0
    cap = arb(".75553")
    assert source[-1].sqrt() < cap
    trace = sum(v.sqrt() for v in source).upper()
    strong = [j for j in range(19) if abs(c[j]) > 4]
    assert strong == [0, 9, 15, 16, 17, 18]
    main = (2, 3, 14)
    others = [j for j in range(19) if j not in strong and j not in main]
    assert len(others) == 10 and all(c[j] < 0 for j in others)
    absmax = max((-c[j]).upper() for j in others)
    errors = lib("fifth-seventh-cutoff-separation")["moment_errors"]
    m2, m4 = errors(
        fields[0], heat, c, [".001" if j in strong else ".75553" for j in range(19)]
    )
    _, other_goal = errors(fields[1], heat, e, [".75553"] * 19)
    slope = -sum(g for g, _, _ in fields[1]) / (8 * arb.pi())
    assert slope > 0
    # 交换第七壳的两条谱线仅置换末端读数。取反后可复用同一个上界证明。
    if swap_five_sixteen:
        e = [-w for w in e]
        e[4], e[15] = e[15], e[4]
    # 由第五壳已证的强层低截止，每个强层对二、四次振幅矩的影响单独包围。
    tv = sum(abs(c[j]) for j in strong)
    low = arb(".001")
    offset2 = m2.upper() + tv * low**2
    offset4 = -m4.lower() + tv * low**4
    rows = []
    for i in range(19):
        h = [abs(overlap[i, j]) for j in range(19)]
        lower, higher = h, h
        if sixth_schur:
            # 保持每个第六壳谱投影的信道可附加任意十九条单位环境向量。
            # 每条路径的复系数固定，环境向量之和的范数受正反三角不等式控制。
            lower, higher = [], []
            for j in range(19):
                lengths = [abs(outgoing[i,k] * incoming[k,j]) for k in range(19)]
                reach = sum(lengths).upper()
                higher.append(min(arb(1), reach))
                lower.append(max(arb(0), (2*max(v.lower() for v in lengths)-reach).lower()))
        elif sixth_measurement:
            # 完全非选择测量的直接和表示给出单个入射层的精确边长平方。
            measured = [
                sum(abs(outgoing[i, k] * incoming[k, j]) ** 2 for k in range(19))
                .nonnegative_part()
                .sqrt()
                for j in range(19)
            ]
            assert all(v.is_finite() for v in measured)
            lower = [a.min(b).lower() for a, b in zip(h, measured, strict=True)]
            higher = [a.max(b).upper() for a, b in zip(h, measured, strict=True)]
        s = sum(higher[j] for j in strong) * low.sqrt()
        a = sum(higher[j] ** 2 for j in others).sqrt()
        rows.append(
            ([lower[j] for j in main], [higher[j] for j in main], s.upper(), a.upper())
        )

    def excluded(left, right):
        tail = trace - sum(left)
        if tail < 0:
            return True
        # 其余十个小负权层的一次和二次平方截止矩 A、B。
        first = -offset2
        second = offset4
        for j, l, h in zip(main, left, right, strict=True):
            first += c[j] * (l if c[j] > 0 else h) ** 2
            second += c[j] * (h if c[j] > 0 else l) ** 4
        if (
            first > absmax * tail**2
            or second < 0
            or (first > 0 and first**3 > second * (absmax * tail) ** 2)
        ):
            return True
        # 对 A 的平方根作用于实际基向量，边长是第五壳对角振幅的平方根。
        rootlo = [v.sqrt() for v in left]
        roothi = [v.sqrt() for v in right]
        root_tail = tail.upper().sqrt()
        score = arb(0)
        for weight, (lower, higher, s, norm) in zip(e, rows, strict=True):
            top = [v * h for v, h in zip(higher, roothi, strict=True)]
            rest = s + norm * root_tail
            if weight > 0:
                y = cap.min((sum(top) + rest) ** 2 + radius)
            else:
                assert weight < 0
                a = max(
                    (v * l - sum(top[j] for j in range(3) if j != i) - rest).lower()
                    for i, (v, l) in enumerate(zip(lower, rootlo, strict=True))
                )
                y = (max(arb(0), a) ** 2 - radius).max(arb(0))
            score += weight * y**4
        return score < arb("-.01")

    floor = tuple(arb(x).sqrt().lower() for x in (".025", ".10", ".001"))
    stack = [(floor, (cap.upper(),) * 3)]
    cells = 0
    while stack:
        left, right = stack.pop()
        cells += 1
        if excluded(left, right):
            continue
        axis = max(range(3), key=lambda i: (right[i] - left[i]).mid())
        assert right[axis] - left[axis] > arb(".0001")
        mid = (left[axis] + right[axis]) / 2
        stack.append((left, right[:axis] + (mid,) + right[axis + 1 :]))
        stack.append((left[:axis] + (mid,) + left[axis + 1 :], right))
    # other_goal = W_target/slope + 对称的热核误差，用其下端抽出有质量 W 的上界。
    if swap_five_sixteen:
        lower = arb("1e-122") + slope * (arb(".01") - other_goal.upper())
        print("交换第 5、16 层后的闭盒数、实际真空项下界", cells, lower.lower().str(15), flush=True)
        assert lower > arb(".09689")
        return lower.lower()
    upper = arb("1e-122") + slope * (arb("-.01") - other_goal.lower())
    print(
        "三维闭盒数、第七壳无质量归一四次矩上界、实际真空 W 上界",
        cells,
        "-.01",
        upper.upper().str(15),
        flush=True,
    )
    assert upper < arb("-.09689")
    if radius > 0:
        gap = source[-1].sqrt() - source[0].sqrt()
        print(
            "允许的逐层读数半径、按原谱跨度换算的酉位移下界",
            radius.str(12),
            (radius / gap).lower().str(14),
            flush=True,
        )
    return upper.upper()


if __name__ == "__main__":
    from argparse import ArgumentParser

    cli = ArgumentParser(description=__doc__)
    cli.add_argument(
        "--sixth-measurement", action="store_true", help="覆盖所有第六壳非选择测量强度"
    )
    cli.add_argument(
        "--readout-radius", default="0", help="允许第七壳各对角读数偏移的半径"
    )
    cli.add_argument("--sixth-schur", action="store_true", help="覆盖保持第六壳各谱层占据数的所有量子信道")
    cli.add_argument("--swap-five-sixteen", action="store_true", help="交换第七壳的第 5、16 层并证明真空项为正")
    options = cli.parse_args()
    verify(
        Path(__file__).resolve().parent,
        sixth_measurement=options.sixth_measurement,
        sixth_schur=options.sixth_schur,
        swap_five_sixteen=options.swap_five_sixteen,
        readout_radius=options.readout_radius,
    )
