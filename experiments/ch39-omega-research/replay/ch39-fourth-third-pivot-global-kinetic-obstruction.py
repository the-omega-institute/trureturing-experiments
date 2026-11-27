"""第四壳第三枢轴五式与全程正动能的六分支连续排除。"""

from pathlib import Path
from runpy import run_path
import inspect

import numpy as np
from flint import arb, ctx

SCALE = 200_000_000


def build(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    number = 2_306_888_100_089
    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=number - 3
    )[2]
    native = lib("native-rh-ray-pvm")
    _, _, pivots, _ = next(
        row for row in native["shell_weights"](root, number) if row[0] == 6
    )
    weights = [v.real for v in pivots[2]]
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root, "1.176397e14", number=number, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, number)
    ctx.prec = 160
    p0, p1 = data["p0"], data["p1"]
    rho = (arb("246.22") / arb("1.22089e19")) ** 2 / 2
    species, _ = neutral["response_rows"](data, masses, signs, rho, arb(".42995401872"))
    energies = [5600, 560000, 56000000, 5600000000, 560000000000]
    pvalues = [(arb(e) / arb("1.22089e19")) ** 2 for e in energies] + [p1]
    # 在低端加密对数截止，五个对数单位后改用 1/32 的步长。
    ticks = np.r_[
        np.arange(125, 5 * 512, 2), np.arange(5 * 512, 80 * 512 + 1, 16)
    ].astype(int)
    assert p0 * (arb(int(ticks[0])) / 512).exp() < ordered[0]
    assert p0 * (arb(int(ticks[-1])) / 512).exp() > ordered[-1]
    high = int(np.searchsorted(ticks, float((arb("1e-6") / p0).log()) * 512) - 1)
    assert p0 * (arb(int(ticks[high])) / 512).exp() < arb("1e-6")
    a1, a2, a3 = data["initial"]
    rows = [
        (
            k,
            z,
            [
                (q - 8 * col / arb(3)) / (4 * arb.pi() * 12),
                (mix - col) / (4 * arb.pi() * 3),
                q / (4 * arb.pi() * 10),
                col / (4 * arb.pi() * 10),
            ],
        )
        for k, z, (q, mix, col) in species
    ]
    initial_readings = [
        (5 * a1 / 3 + a2 - 8 * a3 / 3) / 12,
        (a2 - a3) / 3,
        arb(0),
        *[v for e in energies for v in ((5 * a1 / 3 + a2) / 10, a3 / 10)],
    ]
    kinds = [
        (len(pvalues) - 1, 0),
        (len(pvalues) - 1, 1),
        *[(j, k) for j in range(len(energies)) for k in (2, 3)],
    ]
    norm = [sum(abs(b[k]) for _, _, b in rows) for j, k in kinds]
    rule = gauge["angular_rule"](18, 48)
    weak = native["weak_series"](56)
    # 取原始 56 阶多项式的同一组系数，逐项求其二、三阶导数。
    coefficients = inspect.getclosurevars(weak).nonlocals["coefficients"]
    assert len(coefficients) == 56

    def weak_at(x):
        return weak(x / p0) if x < 1e7 * p0 else neutral["weak_column"](x, p0, p1)

    def point(t):
        return p0 * (arb(int(t)) / 1024).exp()

    start = point(2 * ticks[0])

    def initial():
        raw = [
            [
                gauge["angular_exposure"](k, z, start, p0, p, (p / p0).log(), rule)
                for k, z, b in rows
            ]
            for p in pvalues
        ]
        return [
            -sum(b[k] * v for (_, _, b), v in zip(rows, raw[j], strict=True))
            for j, k in kinds
        ]

    f = initial()
    means = lib("fourth-third-pivot-cyclic-five")["angular_means"]

    def deriv(x):
        sums = [[arb(0) for j in range(4)] for k in range(4)]
        for k, z, b in rows:
            e = (-z / x).exp()
            for j in range(4):
                sums[k][j] += b[j] * e
        base = means(p0 / x)
        diff = [
            [a - b for a, b in zip(base, means(p / x), strict=True)] for p in pvalues
        ]
        return [-sum(sums[k][q] * diff[j][k] for k in range(4)) for j, q in kinds]

    def insert(v, d):
        return [*v[:2], d, *v[2:]]

    left_derivative = deriv(start)
    samples = []
    chord_bounds = []
    samples.append(insert(f, 100 * weak_at(start)))
    for a, b in zip(ticks[:-1], ticks[1:], strict=True):
        x = point(2 * b)
        dl = arb(int(b - a)) / 512
        half = deriv(point(a + b))
        right_derivative = deriv(x)
        # 在对数截止上作辛普森积分。两个能标的指数项各有四阶导数绝对界 17。
        f = [
            v + dl * (l + 4 * m + r) / 6 + arb(0, (dl**5 * 36 * n / 2880).upper())
            for v, l, m, r, n in zip(
                f, left_derivative, half, right_derivative, norm, strict=True
            )
        ]
        # 指数项的二阶导数不超过 4/e²+1/e<1，故响应的三阶导数不超过 2C。
        # 响应二阶导数的均值和离均值余量一起控制整条弦。
        gauge_curvature = [
            abs((r - l) / dl) + 2 * n * dl
            for l, r, n in zip(left_derivative, right_derivative, norm, strict=True)
        ]
        inv = (arb(-int(a)) / 512).exp()
        if inv > arb("1e-7"):
            w2 = abs(sum(n * n * v * inv**n for n, v in enumerate(coefficients, 1)))
            w3 = sum(n**3 * abs(v) * inv**n for n, v in enumerate(coefficients, 1))
            weak_error = (w2 + w3 * dl) * dl**2 / 8 + arb("1e-14")
        else:
            weak_error = 2 * (arb("2.398e-21") + arb("7.468e-36") * inv / p0)
        # 在跨过 1e7 p0 的区间亦可使用级数，其有限高能端误差 2x/p1 < 2e-17。
        chord_errors = insert(
            [v * dl**2 / 8 for v in gauge_curvature], 100 * weak_error
        )
        chord_bounds.append([v.upper() for v in chord_errors])
        left_derivative = right_derivative
        samples.append(insert(f, 100 * weak_at(x)))
    # 每个顶点额外计入左、右两条弦的最大余项。缩放后向外取整。
    lower = []
    upper = []
    padding = [arb(0)] * len(initial_readings)
    for vals, left, right in zip(
        samples, [padding, *chord_bounds], [*chord_bounds, padding], strict=True
    ):
        bounded = [
            v + arb(0, max(a, b).upper())
            for v, a, b in zip(vals, left, right, strict=True)
        ]
        lower.append([int((v * SCALE).lower().floor().unique_fmpz()) for v in bounded])
        upper.append([int((v * SCALE).upper().ceil().unique_fmpz()) for v in bounded])
    return (
        np.array(lower, np.int64),
        np.array(upper, np.int64),
        high,
        weights,
        initial_readings,
    )


def verify(root):
    gravity = run_path(str(root / "ch39-fourth-third-pivot-strong-gravity-cutoff.py"))
    gravity["verify"](root)
    lower, upper, high, weights, start = build(root)
    strong = [i for i, c in enumerate(weights) if abs(c) > 4]
    assert [i + 1 for i in strong] == [1, 10, 16, 17, 18, 19]
    source = root / "ch39-fourth-third-pivot-global-kinetic-duals.txt"
    lines = iter(source.read_text().splitlines())
    results = []
    bound = max(int(np.max(np.abs(lower))), int(np.max(np.abs(upper))))
    for distinguished in strong:
        assert next(lines) == "branch " + str(distinguished + 1)
        full = [(0, len(lower) - 1) for _ in range(19)]
        full[distinguished] = (high, len(lower) - 1)
        stack = [full]
        leaves = 0
        margin = arb(100)
        while stack:
            box = stack.pop()
            operation, *arguments = next(lines).split()
            if operation == "split":
                axis, mid = map(int, arguments)
                axis -= 1
                assert 0 <= axis < 19
                lo, hi = box[axis]
                assert lo < mid < hi
                stack.append([*box[:axis], (mid, hi), *box[axis + 1 :]])
                stack.append([*box[:axis], (lo, mid), *box[axis + 1 :]])
                continue
            assert operation == "dual"
            multipliers = [0] * 16
            for argument in arguments:
                i, value = map(int, argument.split("="))
                assert 0 <= i < 16 and value >= 0 and multipliers[i] == 0
                multipliers[i] = value
            # H=(R,-R,-A)。相同 R 的两个方向相减，所有动力学行以非正号出现。
            direction = [multipliers[i] - multipliers[i + 3] for i in range(3)] + [
                -v for v in multipliers[6:]
            ]
            assert bound * sum(map(abs, direction)) < 2**62
            pos = np.array([max(v, 0) for v in direction], np.int64)
            neg = np.array([min(v, 0) for v in direction], np.int64)
            floor = lower @ pos + upper @ neg
            ceil = upper @ pos + lower @ neg
            total = sum(a * d for a, d in zip(start, direction, strict=True)) * SCALE
            for c, (lo, hi) in zip(weights, box, strict=True):
                minimum = int(np.min(floor[lo : hi + 1]))
                maximum = int(np.max(ceil[lo : hi + 1]))
                # 先按权的严格符号选取区间端点，再用 Arb 乘实际谱权。
                endpoint = minimum if c > 0 else maximum
                assert c > 0 or c < 0
                total += c * endpoint
            certified = total.lower() / (SCALE * 10**9)
            assert certified > 0, (distinguished + 1, box, certified)
            margin = min(margin, certified.lower())
            leaves += 1
        results.append((leaves, margin))
        print(
            "高截止层、分离叶数、严格最小余量",
            distinguished + 1,
            leaves,
            margin.str(13),
            flush=True,
        )
    assert next(lines, None) is None
    print("六分支合计分离叶数", sum(n for n, _ in results), flush=True)
    return results


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
