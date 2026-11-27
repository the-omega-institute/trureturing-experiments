"""固定第五壳相干解，在第七壳同时恢复引力双矩和弱差所需的最少二层酉门数。"""

from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


def inputs(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    base = lib("fifth-seventh-two-pair-thermal-rotation")
    _, _, _, _, boxes, _ = lib("fifth-third-pivot-three-high-positive")["endpoints"](
        root
    )
    frame = base["coherent_frame"](root, lib, boxes)
    diag, cap = base["readings"](frame)
    number = base["SEVENTH"]
    pivots = next(
        row
        for row in lib("native-rh-ray-pvm")["shell_weights"](root, number)
        if row[0] == 6
    )[2]
    weights = [w.real for w in pivots[2]]
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](
        root, [number]
    )
    ctx.prec = 192
    slope, column = base["thermal_columns"](field, moments, heat)
    fixed = [[w * v for v in column(y)] for w, y in zip(weights, diag, strict=True)]
    return diag, cap, weights, field, heat, fixed, column, slope, frame


def verify(root):
    data = inputs(root)
    diag, cap, weights, field, heat, fixed, column, _, _ = data
    delta = arb("1e-14")
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    correction = sum(map(abs, weights)) * (
        arb("2.398e-21") + arb("7.468e-36") / delta**2
    )
    assert delta**2 > (arb("91.1876") / arb("1.22089e19")) ** 2
    assert all(y > delta for y in diag) and reference + correction < arb("-.074")
    # 若三式同时成立，有一个活动振幅必在 [0,delta]。
    small = arb(0).union(delta)
    x = cap.upper() ** 2
    pi = arb.pi()
    ks = sum(a for _, a, _ in field) / (2 * pi)
    ws = -sum(g for g, _, _ in field) / (8 * pi)
    errors = [
        x * sum(abs(a) * (1 - heat(2, z / x)) for _, a, z in field) / (2 * pi),
        x**2
        * sum(abs(g) * (arb(1) / 2 - heat(3, z / x)) for g, _, z in field)
        / (4 * pi),
    ]
    slack = arb(".003")
    all_layers = set(range(19))
    split_counts = [0, 0]
    supports = [0, 0]
    rigorous = 0

    def nonnegative(lo, hi):
        return max(arb(0), lo.lower()), min(cap.upper(), hi.upper())

    def check(low, rest_pair, partner):
        nonlocal rigorous
        i, j = rest_pair
        topology = int(partner is not None)
        changed = {low, i, j}
        total = diag[i] + diag[j]
        pinned = [(low, small)]
        if partner is None:
            total += diag[low] - small
        else:
            changed.add(partner)
            pinned.append((partner, diag[low] + diag[partner] - small))
        other = all_layers - changed
        targets = (1, arb("1e-122"))
        offset = [sum(fixed[k][n] for k in other) - v for n, v in enumerate(targets)]
        cost = sum(abs(weights[k]) for k in changed)
        goal = [
            -v + arb(0, (slack + cost * e).upper())
            for v, e in zip(offset, errors, strict=True)
        ]
        for n, (scale, p) in enumerate(((ks, 2), (ws, 4))):
            goal[n] -= scale * sum(weights[k] * y**p for k, y in pinned)
        supports[topology] += 1
        stack = [nonnegative(arb(0), total)]
        while stack:
            left, right = stack.pop()
            bl, br = nonnegative(total - right, total - left)
            if bl > br:
                continue
            for n, (scale, p) in enumerate(((ks, 2), (ws, 4))):
                lead = scale * (
                    weights[i] * (left**p).union(right**p)
                    + weights[j] * (bl**p).union(br**p)
                )
                if not lead.overlaps(goal[n]):
                    break
            else:
                if right - left < arb(".0001"):
                    # 在无质量余项过宽的盒内，改用实际有质量核的均值包围。
                    actual = offset.copy()
                    for k, y in [*pinned, (i, left.union(right)), (j, bl.union(br))]:
                        if y.lower() > 0:
                            val = column(y)
                        else:
                            high = y.upper()
                            val = [
                                arb(
                                    0,
                                    (
                                        sum(abs(a) for _, a, _ in field)
                                        * high**2
                                        / (2 * pi)
                                    ).upper(),
                                ),
                                arb(
                                    0,
                                    (
                                        sum(abs(g) for g, _, _ in field)
                                        * high**4
                                        / (8 * pi)
                                    ).upper(),
                                ),
                            ]
                        for n in range(2):
                            actual[n] += weights[k] * val[n]
                    rigorous += 1
                    if any(abs(v) > slack for v in actual):
                        split_counts[topology] += 1
                        continue
                assert right - left > arb("1e-9"), (
                    low,
                    rest_pair,
                    partner,
                    left,
                    right,
                )
                mid = (left + right) / 2
                stack.extend(((left, mid), (mid, right)))
                continue
            split_counts[topology] += 1

    for low in range(19):
        for pair in combinations(sorted(all_layers - {low}), 2):
            check(low, pair, None)
            for partner in sorted(all_layers - {low, *pair}):
                check(low, pair, partner)
    assert supports == [2907, 46512]
    print(
        "相交与不交拓扑的带低层选择数、两类分离闭段数、实际核细化数、引力残差界",
        supports,
        split_counts,
        rigorous,
        slack.str(3),
        flush=True,
    )

    three_gates(root, data)


def three_gates(root, data):
    diag, _, weights, _, _, fixed, column, slope, frame = data
    floor = frame[0]
    weak_faces = (arb("6.5e-17"), arb("9e-17"))
    small = weak_faces[0].union(weak_faces[1])
    tri = diag[4] + diag[9] + diag[15]
    pair = diag[6] + diag[17]
    excluded = {4, 9, 15, 6, 17}
    fixed_sum = [
        sum(row[n] for k, row in enumerate(fixed) if k not in excluded) - target
        for n, target in enumerate((1, arb("1e-122")))
    ]

    def readings(center):
        a, b = center
        return a, tri - a - small, small, b, pair - b

    idx = (4, 9, 15, 6, 17)

    def fun(center):
        vals = [column(y) for y in readings(center)]
        return arb_mat(
            [
                [
                    fixed_sum[n]
                    + sum(weights[k] * v[n] for k, v in zip(idx, vals, strict=True))
                ]
                for n in range(2)
            ]
        )

    def jac(center):
        a, h, _, b, v = readings(center)
        cols = []
        for i, j, y, z in ((4, 9, a, h), (6, 17, b, v)):
            da, db = slope(y), slope(z)
            cols.append([weights[i] * da[n] - weights[j] * db[n] for n in range(2)])
        return arb_mat(cols).transpose()

    center = [arb(".1084934024"), arb(".05367806836")]
    for _ in range(2):
        shift = jac(center).mid().solve(fun(center).mid()).entries()
        center = [(a - d).mid() for a, d in zip(center, shift, strict=True)]
    rad = arb("1e-5")
    box = [v + arb(0, rad.upper()) for v in center]
    pre = jac(center).mid().inv().mid()
    assert not pre.det().contains(0)
    delta = (pre * fun(center)).entries()
    err = (arb_mat([[1, 0], [0, 1]]) - pre * jac(box)).tolist()
    ratio = max(
        (abs(v) / rad + sum(map(abs, row))).upper()
        for v, row in zip(delta, err, strict=True)
    )
    print(
        "三门的第五、第七层振幅中心、半径、预条件面界",
        *[a.str(15) for a in (*center, rad, ratio)],
        flush=True,
    )
    assert ratio < arb(".050005")
    # 移动两个引力目标时，预条件化残差的每个分量至多增加此行的绝对和乘目标半径。
    inverse_norm = max(sum(abs(v) for v in row).upper() for row in pre.tolist())
    target_width = rad*(1-ratio)/inverse_norm
    print("外部两个引力目标可同时变动的半宽下界",target_width.lower().str(15),flush=True)
    assert target_width > arb(".00014")
    rows = readings(box)
    restored = dict(zip(idx, rows, strict=True))
    assert diag[17] < rows[3] < diag[6]
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0 = (arb("91.1876") / arb("1.22089e19")) ** 2
    assert all(restored.get(i, diag[i]) ** 2 > p0 for i in range(19))
    rest_error = sum(
        abs(weights[i])
        * (arb("2.398e-21") + arb("7.468e-36") / restored.get(i, diag[i]) ** 2)
        for i in range(19)
        if i != 15
    )
    weak = run_path(str(root / "ch39-native-rh-ray-pvm.py"))["weak_series"](56)
    left, right = [
        reference
        + weights[15] * (weak(z * z / p0) - reference)
        + arb(0, rest_error.upper())
        for z in weak_faces
    ]
    print("第十六层振幅上下两面的弱差", left.str(12), right.str(12), flush=True)
    assert left > arb(".03408") and right < arb("-.01846")
    # 判断整个振幅盒能否用 (5,10)、(10,16) 两门同时得到第一组三个读数。
    m, h, k, out, overlap, probs, phases = frame

    def minor(i, j):
        return ((h - m) * (k - m)).sqrt() * sum(
            (overlap[i, a] * overlap[j, b] - overlap[i, b] * overlap[j, a])
            * (probs[a] * probs[b]).sqrt()
            * (phases[b] - phases[a])
            for a, b in combinations(range(19), 2)
            if a in (3, 14) or b in (3, 14)
        )

    def adj(matrix):
        return matrix.conjugate().transpose()

    inner = (h - m) * out[4, 0] * out[9, 0].conjugate() + (k - m) * out[4, 1] * out[
        9, 1
    ].conjugate()
    top = (
        diag[4]
        + diag[9]
        - 2 * m
        + ((diag[4] - diag[9]) ** 2 + 4 * abs(inner) ** 2).sqrt()
    ) / 2
    bottom = abs(minor(4, 9)) ** 2 / top
    retained = diag[4] + diag[9] - rows[0] - m
    prob = (retained - bottom) / (top - bottom)
    assert 0 < prob < 1
    proj = acb_mat(
        [[diag[4] - m - bottom, inner], [inner.conjugate(), diag[9] - m - bottom]]
    ) / (top - bottom)
    vec = acb_mat([[minor(4, 15)], [minor(9, 15)]])
    plus = (adj(vec) * proj * vec)[0, 0].real
    minus = (adj(vec) * (acb_mat([[1, 0], [0, 1]]) - proj) * vec)[0, 0].real
    assert plus > 0 and minus > 0
    least = ((prob * plus).sqrt() - ((1 - prob) * minus).sqrt()) ** 2
    low = small - floor
    trace = retained + diag[15] - m
    assert low > 0 and trace > 2 * low
    test = least / (low * (trace - low))
    print("两次相交门的行列式与可达阈值之比上界", test.upper().str(14), flush=True)
    assert test < arb(".009095")
    return [restored.get(i,diag[i]) for i in range(19)]


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
