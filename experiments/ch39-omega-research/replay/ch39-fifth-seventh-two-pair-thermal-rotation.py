"""第五壳相干正解在第七壳以两次不交二层酉旋转恢复双引力矩。"""

from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

FIFTH = 2_306_888_100_090
SEVENTH = 2_306_888_100_092
GAUGE_ROWS = (
    (11, 0, 1, 0, 1, 0, 1, 11, 11, 11, 10, 10, 9, 9, 9, 14, 6, 14, 6),
    (11, 0, 1, 1, 0, 0, 1, 11, 11, 11, 10, 10, 9, 9, 9, 14, 6, 14, 6),
)
PAIRS = ((4, 15), (14, 17))


def adj(matrix):
    return matrix.conjugate().transpose()


def coherent_reading(root, lib, boxes):
    _, source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2]
    low, h, k = [source[i].sqrt() for i in (0, -2, -1)]
    assert 0 < low < h < k and all(x is source[0] for x in source[:17])
    probs = [(x.sqrt() - low) / (h + k - 2 * low) for x in boxes]
    assert all(p > 0 for p in probs) and sum(probs).contains(1)
    a, b = 3, 14
    rest = 1 - probs[a] - probs[b]
    cosine = (rest**2 - probs[a] ** 2 - probs[b] ** 2) / (2 * probs[a] * probs[b])
    assert rest > 0 and -1 < cosine < 1
    v = acb(cosine, (1 - cosine**2).sqrt())
    phases = [-(probs[a] + probs[b] * v) / rest] * 19
    phases[a], phases[b] = acb(1), v
    cols = acb_mat(
        [[p.sqrt(), p.sqrt() * phase] for p, phase in zip(probs, phases, strict=True)]
    )
    gram = adj(cols) * cols
    assert all(gram[i, j].contains(int(i == j)) for i in range(2) for j in range(2))
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    bases = []
    for number, rows in zip((FIFTH, SEVENTH), GAUGE_ROWS, strict=True):
        basis = spectrum(shell, 6, number)[1]
        for j, i in enumerate(rows):
            entry = basis[i, j]
            assert abs(entry) > 0
            phase = entry.conjugate() / abs(entry)
            for m in range(19):
                basis[m, j] *= phase
        bases.append(basis)
    out = adj(bases[1]) * bases[0] * cols
    diag = [
        low + (h - low) * abs(out[i, 0]) ** 2 + (k - low) * abs(out[i, 1]) ** 2
        for i in range(19)
    ]
    assert all(y > 0 for y in diag) and sum(diag).overlaps(17 * low + h + k)
    return diag, k


def thermal_columns(field, moments, heat):
    def slope(y):
        exps = [(-z / y**2).exp() for _, _, z in field]
        return (
            y * sum(a * e for (_, a, _), e in zip(field, exps, strict=True)) / arb.pi(),
            -(y**3)
            * sum(g * e for (g, _, _), e in zip(field, exps, strict=True))
            / (2 * arb.pi()),
        )

    def column(y):
        mid = y.mid()
        base = moments(field, heat, mid**2)
        return [v + (y - mid) * d for v, d in zip(base, slope(y), strict=True)]

    return slope, column


def two_rotations(diag, weights, fixed, slope, column):
    active = {i for pair in PAIRS for i in pair}
    rest = [
        sum(row[n] for i, row in enumerate(fixed) if i not in active) - target
        for n, target in enumerate((1, arb("1e-122")))
    ]

    def readings(i, j, t):
        diff = diag[j] - diag[i]
        return diag[i] + t * diff, diag[j] - t * diff

    def residual(params):
        sums = rest.copy()
        for (i, j), t in zip(PAIRS, params, strict=True):
            vi, vj = [column(y) for y in readings(i, j, t)]
            for n in range(2):
                sums[n] += weights[i] * vi[n] + weights[j] * vj[n]
        return arb_mat([[v] for v in sums])

    def jac(params):
        cols = []
        for (i, j), t in zip(PAIRS, params, strict=True):
            vi, vj = [slope(y) for y in readings(i, j, t)]
            cols.append(
                [
                    (diag[j] - diag[i]) * (weights[i] * vi[n] - weights[j] * vj[n])
                    for n in range(2)
                ]
            )
        return arb_mat(cols).transpose()

    center = [arb(".42757613"), arb(".40844060")]
    for _ in range(2):
        shift = jac(center).mid().solve(residual(center).mid()).entries()
        center = [(t - d).mid() for t, d in zip(center, shift, strict=True)]
    radius = arb("1e-4")
    box = [t + arb(0, radius.upper()) for t in center]
    assert all(0 < t < 1 for t in box)
    pre = jac(center).mid().inv().mid()
    assert not pre.det().contains(0)
    shift = (pre * residual(center)).entries()
    error = arb_mat([[1, 0], [0, 1]]) - pre * jac(box)
    ratio = max(
        (abs(v) / radius + sum(map(abs, row))).upper()
        for v, row in zip(shift, error.tolist(), strict=True)
    )
    print(
        "两个旋转概率中心、共同半径、二维自映比",
        *(v.str(14) for v in center),
        radius.str(4),
        ratio.str(14),
        flush=True,
    )
    assert ratio < 1
    return box


def single_rotation(diag, weights, fixed, cap, field, heat):
    # 正核的质量误差随平方截止单调，由源谱上界封住两条改变量。
    x = cap.upper() ** 2
    pi = arb.pi()
    slopes = (
        sum(a for _, a, _ in field) / (2 * pi),
        -sum(g for g, _, _ in field) / (8 * pi),
    )
    errors = (
        x * sum(abs(a) * (1 - heat(2, z / x)) for _, a, z in field) / (2 * pi),
        x**2
        * sum(abs(g) * (arb(1) / 2 - heat(3, z / x)) for g, _, z in field)
        / (4 * pi),
    )
    slack = arb(".35")
    npairs, leaves = 0, 0
    for i, j in combinations(range(19), 2):
        npairs += 1
        total = diag[i] + diag[j]
        coeff = abs(weights[i]) + abs(weights[j])
        rest = [
            sum(row[n] for k, row in enumerate(fixed) if k not in (i, j)) - target
            for n, target in enumerate((1, arb("1e-122")))
        ]
        goals = [
            -r + arb(0, (slack + coeff * e).upper())
            for r, e in zip(rest, errors, strict=True)
        ]
        stack = [(arb(0), min(cap.upper(), total.upper()))]
        while stack:
            l, r = stack.pop()
            bl = max(arb(0), (total - r).lower())
            br = min(cap.upper(), (total - l).upper())
            if bl > br:
                continue
            for n, power in enumerate((2, 4)):
                approx = slopes[n] * (
                    weights[i] * (l**power).union(r**power)
                    + weights[j] * (bl**power).union(br**power)
                )
                if not approx.overlaps(goals[n]):
                    leaves += 1
                    break
            else:
                assert r - l > arb("1e-9"), (i + 1, j + 1, l, r, goals)
                m = (l + r) / 2
                stack.extend(((l, m), (m, r)))
    assert npairs == 171
    print(
        "任意单旋转的层对数、分离盒数、双矩残差最大值严格下界",
        npairs,
        leaves,
        slack.str(3),
        flush=True,
    )


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    _, _, _, _, boxes, _ = lib("fifth-third-pivot-three-high-positive")["endpoints"](
        root
    )
    diag, cap = coherent_reading(root, lib, boxes)
    pivots = next(
        row
        for row in lib("native-rh-ray-pvm")["shell_weights"](root, SEVENTH)
        if row[0] == 6
    )[2]
    weights = [w.real for w in pivots[2]]
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](
        root, [SEVENTH]
    )
    ctx.prec = 192
    slope, column = thermal_columns(field, moments, heat)
    fixed = [[w * v for v in column(y)] for w, y in zip(weights, diag, strict=True)]
    two_rotations(diag, weights, fixed, slope, column)
    single_rotation(diag, weights, fixed, cap, field, heat)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
