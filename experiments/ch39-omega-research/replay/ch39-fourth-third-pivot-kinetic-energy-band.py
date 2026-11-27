"""改变匹配能标并逐点重选截止时，第四壳第三枢轴的正动能排除窗。"""

from pathlib import Path
from runpy import run_path
import json
import numpy as np
from flint import arb


def refined_boxes(box, original, operations, decode):
    if operations is None:
        yield box, original
        return
    stack = [box]
    for action, *values in operations:
        section = stack.pop()
        if action == "s":
            axis, mid = values
            assert 0 <= axis < 19
            lo, hi = section[axis]
            assert lo < mid < hi
            stack.extend(
                (
                    [*section[:axis], (mid, hi), *section[axis + 1 :]],
                    [*section[:axis], (lo, mid), *section[axis + 1 :]],
                )
            )
        elif action == "o":
            yield section, original
        else:
            assert action == "d"
            multipliers = [0] * 16
            for index, integer in values[0]:
                assert 0 <= index < 16 and integer >= 0 and multipliers[index] == 0
                multipliers[index] = integer
            yield section, decode(multipliers)
    assert not stack


def second_pivot(root, core, order, tables, minimum):
    lower, upper, _, _, initial, coeff = tables
    verify = run_path(str(root / "ch39-fourth-second-pivot-obstruction.py"))["verify"]
    weights, threshold, cut, p0, weak_cap = verify(root, weak_minimum=minimum)
    ticks = core["cutoff_ticks"]()
    low_end = int(np.searchsorted(ticks, float((threshold / p0).log() * 512))) + 1
    high_start = int(np.searchsorted(ticks, float((cut / p0).log() * 512))) - 1
    assert p0 * (arb(int(ticks[low_end])) / 512).exp() > threshold
    assert p0 * (arb(int(ticks[high_start])) / 512).exp() < cut
    variation = sum(map(abs, weights))
    width = max(
        b - a for a, b in zip(sorted(order)[:-1], sorted(order)[1:], strict=True)
    )
    chord = (
        width**2
        * variation
        * sum(abs(a) for a, _ in coeff)
        / (arb(1).exp() * 8 * 64**2)
    )
    bounds = []
    for small in (0, 15):
        box = [(0, len(lower) - 1) for _ in weights]
        box[small] = (0, low_end)
        box[17] = (high_start, len(lower) - 1)
        projected = []
        for j in range(len(order)):
            col = 0 if j == 0 else 13 + 2 * (j - 1)
            projected.append(
                core["signed_lower"](
                    lower[:, col],
                    upper[:, col],
                    box,
                    weights,
                    initial[0] * core["SCALE"],
                )
                / core["SCALE"]
            )
        bounds.append(12 * (min(projected) - chord))
    assert all(v > 0 for v in bounds)
    print(
        "第二枢轴弱差上界及两个低层分支在整窗的光子残差下界",
        weak_cap.upper().str(13),
        *(v.lower().str(13) for v in bounds),
        flush=True,
    )
    return bounds


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    lib("fourth-third-pivot-strong-gravity-cutoff")["verify"](root)
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    specification = json.loads(
        (root / "ch39-fourth-third-pivot-kinetic-energy-band-duals.json").read_text()
    )
    windows = specification["windows"]
    lo, hi = windows[0]["low"], windows[-1]["high"]
    assert lo <= 0 <= hi
    assert all(
        a["high"] == b["low"] for a, b in zip(windows[:-1], windows[1:], strict=True)
    )
    for window in windows:
        assert window["low"] < window["high"] and window["step"] > 0
        assert (window["high"] - window["low"]) % window["step"] == 0
    nodes = sorted(
        {i for w in windows for i in range(w["low"], w["high"] + 1, w["step"])}
    )
    order = [0, *[i for i in nodes if i != 0]]
    tables = core["build_grid"](root, order)
    lower, upper, high, weights, start, coeff = tables
    edges = [arb("1.176397e14") * (arb(j) / 128).exp() for j in (lo, hi)]
    # 第一枢轴有独立的连续能窗 [1e13,1e17] 证书，第二枢轴在此窗逐层重算。
    first = lib("fourth-variable-scale-neutral-gauge")["verify"](root)
    assert arb("1e13") < edges[0] < edges[1] < arb("1e17") and all(
        r[0] > 0 for r in first
    )
    second_pivot(root, core, order, tables, edges[0].lower())
    scale = core["SCALE"]
    bound = max(int(np.max(abs(lower))), int(np.max(abs(upper))))
    variation = sum(map(abs, weights))
    results = {}
    for window in windows:
        selected = [
            order.index(k)
            for k in range(window["low"], window["high"] + 1, window["step"])
        ]
        replacements = dict(window["leaves"])
        seen = set()
        for i, (axis, box, original) in enumerate(
            core["source_boxes"](root, high, len(lower) - 1, weights)
        ):
            if i in replacements:
                seen.add(i)
            for section, direction in refined_boxes(
                box, original, replacements.get(i), core["decode_direction"]
            ):
                assert bound * sum(map(abs, direction)) < 2**62
                # 固定的五处动力学读数和全窗一致的弱差盒，只需作一次顶点点积。
                floor, ceil = core["projection"](
                    lower[:, 2:13], upper[:, 2:13], direction[2:], bound
                )
                initial = (
                    sum(a * d for a, d in zip(start, direction, strict=True)) * scale
                )
                row = []
                for j in selected:
                    start_column = 0 if j == 0 else 13 + 2 * (j - 1)
                    columns = slice(start_column, start_column + 2)
                    pf, pc = core["projection"](
                        lower[:, columns], upper[:, columns], direction[:2], bound
                    )
                    row.append(
                        core["signed_lower"](
                            floor + pf, ceil + pc, section, weights, initial
                        )
                        / (scale * 10**9)
                    )
                # 沿 log p 插值。每个物种的二阶导数绝对值至多 1/e。
                projected = sum(
                    abs(direction[0] * a + direction[1] * b) for a, b in coeff
                )
                chord = (
                    window["step"] ** 2
                    * variation
                    * projected
                    / (arb(1).exp() * 8 * 64**2 * 10**9)
                )
                margin = (min(row) - chord).lower()
                assert margin > 0, (i, axis + 1, margin)
                count, old = results.get(axis, (0, arb(100)))
                results[axis] = count + 1, min(margin, old)
        assert set(replacements) == seen
    for axis, (count, margin) in results.items():
        print(
            "高层、闭盒数、整个匹配能窗的最小严格余量",
            axis + 1,
            count,
            margin.str(11),
            flush=True,
        )
    ends = [arb("1.176397e14") * (arb(j) / 128).exp() for j in (lo, hi)]
    print(
        "能窗两端、全部闭盒数",
        *[v.str(13) for v in ends],
        sum(n for n, _ in results.values()),
        flush=True,
    )
    return ends, results


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
