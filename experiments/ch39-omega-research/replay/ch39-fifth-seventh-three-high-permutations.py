"""第五壳最小三高解的全部十九层重排，与第七壳的引力双矩严格分离。"""

from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def thermal_bounds(field, heat, cap):
    """归一到无质量斜率后的大位置质量误差与小位置绝对读数上界。"""
    pi = arb.pi()
    ka = sum(a for _, a, _ in field) / (2 * pi)
    wa = -sum(g for g, _, _ in field) / (8 * pi)
    assert ka < 0 < wa
    high = (
        cap
        * sum(abs(a) * (1 - heat(2, z / cap)) for _, a, z in field)
        / (2 * pi * abs(ka)),
        cap**2
        * sum(abs(g) * (arb(1) / 2 - heat(3, z / cap)) for g, _, z in field)
        / (4 * pi * wa),
    )
    # 正核 x E₂(z/x) 与 x² E₃(z/x) 单调。保留小位置各场的抵消。
    previous = [(arb(0), arb(0)) for _ in field]
    low = [arb(0), arb(0)]
    for j in range(1, 65):
        x = arb("1e-6") * j / 64
        current = [(x * heat(2, z / x), x * x * heat(3, z / x)) for _, _, z in field]
        k = sum(
            a * left[0].union(right[0])
            for (_, a, _), left, right in zip(field, previous, current, strict=True)
        ) / (2 * pi * ka)
        w = -sum(
            g * left[1].union(right[1])
            for (g, _, _), left, right in zip(field, previous, current, strict=True)
        ) / (4 * pi * wa)
        low = [
            max(old, abs(value).upper()) for old, value in zip(low, (k, w), strict=True)
        ]
        previous = current
    return ka, wa, high, low


def permitted(weights, indices, bounds, slack):
    ka, wa, high, low = bounds
    chosen = [weights[i] for i in indices]
    large = sum(map(abs, chosen))
    small = sum(abs(w) for i, w in enumerate(weights) if i not in indices)
    goals = []
    for anchor, scale, eh, el in zip(
        (1, arb("1e-122")), (ka, wa), high, low, strict=True
    ):
        err = large * eh + small * el + slack / abs(scale)
        goals.append(anchor / scale + arb(0, err.upper()))
    return chosen, goals


def intersects(candidate, second, fourth):
    weights, goals = candidate
    for powers, goal in zip((second, fourth), goals, strict=True):
        if not sum(w * y for w, y in zip(weights, powers, strict=True)).overlaps(goal):
            return False
    return True


def verify(root, *, slack=".04"):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    numbers = (2_306_888_100_090, 2_306_888_100_092)
    _, cutoffs = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2]
    native = lib("native-rh-ray-pvm")
    weights = []
    for number in numbers:
        pivots = next(
            row for row in native["shell_weights"](root, number) if row[0] == 6
        )[2]
        weights.append([w.real for w in pivots[2]])
    _, heat, fields = lib("three-index-free-atoms")["broken_fields"](root, numbers)
    ctx.prec = 192
    trace = sum(x.sqrt() for x in cutoffs)
    ymax = cutoffs[-1].sqrt().upper()
    trace_min = (trace - 16 * arb(".001")).lower()
    trace_max = trace.upper()
    bounds = [thermal_bounds(rows, heat, ymax**2) for rows in fields]
    fifth = permitted(weights[0], (2, 3, 14), bounds[0], arb(0))
    options = [
        permitted(weights[1], indices, bounds[1], arb(slack))
        for indices in permutations(range(19), 3)
    ]
    assert len(options) == 5814
    # 覆盖整个振幅立方，而不预选热矩方程的任何一支根。
    stack = [((arb(0),) * 3, (ymax,) * 3, options)]
    nodes = 0
    while stack:
        lower, upper, candidates = stack.pop()
        nodes += 1
        lower = tuple(
            max(v, (trace_min - sum(upper[j] for j in range(3) if j != i)).lower())
            for i, v in enumerate(lower)
        )
        upper = tuple(
            min(v, (trace_max - sum(lower[j] for j in range(3) if j != i)).upper())
            for i, v in enumerate(upper)
        )
        if any(a > b for a, b in zip(lower, upper, strict=True)):
            continue
        second = [(a * a).union(b * b) for a, b in zip(lower, upper, strict=True)]
        fourth = [(a**4).union(b**4) for a, b in zip(lower, upper, strict=True)]
        if not intersects(fifth, second, fourth):
            continue
        candidates = [one for one in candidates if intersects(one, second, fourth)]
        if not candidates:
            continue
        axis = max(range(3), key=lambda i: (upper[i] - lower[i]).mid())
        assert upper[axis] - lower[axis] > arb("1e-8")
        mid = (lower[axis] + upper[axis]) / 2
        stack.extend(
            (
                (lower, upper[:axis] + (mid,) + upper[axis + 1 :], candidates),
                (lower[:axis] + (mid,) + lower[axis + 1 :], upper, candidates),
            )
        )
    print(
        "高层有序排列数、闭盒数、第七壳两引力残差最大值严格下界",
        len(options),
        nodes,
        slack,
        flush=True,
    )
    return nodes


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
