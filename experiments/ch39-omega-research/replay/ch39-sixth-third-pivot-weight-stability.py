"""第六壳第三枢轴的五式与全程正性在归一化谱权扰动下仍不相容。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb


NUMBER = 2_306_888_100_091
BRANCHES = (0, 9, 15)
DUALS = "ch39-sixth-third-pivot-global-kinetic-duals.txt"


def least_sum(weights, bounds, radius):
    """各层响应落在给定区间时，对归一权向量的整个闭球取下界。"""
    half = radius / 2

    def term(weight, lo, hi):
        return (weight * lo).min(weight * hi)

    center = []
    gains = []
    losses = []
    for weight, (lo, hi) in zip(weights, bounds, strict=True):
        value = term(weight, lo, hi)
        center.append(value)
        gains.append(term(weight + half, lo, hi) - value)
        losses.append(term(weight - half, lo, hi) - value)
    # 每个顶点都从一层移走半径的一半，再把它加至另一层。
    change = min((gain + loss).lower() for i, gain in enumerate(gains)
                 for j, loss in enumerate(losses) if i != j)
    return sum(center).lower() + min(arb(0), change)


def gravity_bounds(root):
    atoms = run_path(str(root / "ch39-three-index-free-atoms.py"))
    source = run_path(str(root / "ch39-adjacent-ordered-charge-six-cutoffs.py"))
    _, heat, (field,) = atoms["broken_fields"](root, [NUMBER])
    _, ordered = source["verify"](root, refine_node=2_306_888_100_086)[2]
    threshold, cutoff = arb("1e-6"), ordered[-1].upper()
    assert threshold < cutoff < 1
    grid = [arb(2)**(arb(k)/64) for k in range(-1920, 1)]
    edges = sorted([arb(0), threshold, cutoff, *(v for v in grid if v < cutoff)],
                   key=lambda v: v.mid())
    pi = arb.pi()

    def energy(pairs):
        return sum(a*k/(2*pi) - 250*g*w/(4*pi)
                   for (g, a, _), (k, w) in zip(field, pairs, strict=True))

    low = []
    all_values = []
    last_value = arb(0)
    for left, right in zip(edges, edges[1:]):
        assert left < right
        kernels = [(right*heat(2, z/right), right**2*heat(3, z/right)) for _, _, z in field]
        value = energy(kernels)
        if left > 0:
            cell = left.union(right)
            slope = sum((a/(2*pi) - 250*g*cell/(4*pi)) * (-z/cell).exp() for g, a, z in field)
            shift = arb(0).union(right-left) * slope
            band = (last_value + shift).intersection(value - shift)
        else:
            # 两个热核都从零单调增加，首格包含原点。
            band = energy([(k.union(arb(0)), w.union(arb(0))) for k, w in kernels])
        all_values.append(band)
        if left < threshold:
            low.append(band)
        last_value = value

    def negative_range(values):
        return -max(x.upper() for x in values), -min(x.lower() for x in values)

    # 用 -F 的下界得到 F 的上界。
    full_band, low_band = negative_range(all_values), negative_range(low)
    return [low_band if i in BRANCHES else full_band for i in range(19)]


def verify(root):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    lower, upper, high, weights, start = core["build_grid"](root, (0,), number=NUMBER)[:5]
    assert len(weights) == 19 and sum(weights).contains(1)
    radius = arb(".030568")
    ceiling = -least_sum(weights, gravity_bounds(root), radius)
    print("三个指定层全低时 K+250W 的上界", ceiling.upper().str(14), flush=True)
    assert ceiling < arb(".9999")

    max_entry = max(int(np.max(np.abs(lower))), int(np.max(np.abs(upper))))
    scale = core["SCALE"]
    assert scale % 10 == 0
    branch_radius = dict(zip(BRANCHES, map(arb, (".36", ".12", ".16")), strict=True))
    counts = dict.fromkeys(BRANCHES, 0)
    deficits = {axis: [] for axis in BRANCHES}
    for axis, box, direction in core["source_boxes"](root, high, len(lower)-1, weights, filename=DUALS):
        if axis not in branch_radius:
            continue
        floor, ceil = core["projection"](lower, upper, direction, max_entry)
        bands = [(int(floor[a:b+1].min()), int(ceil[a:b+1].max())) for a, b in box]
        constant = scale * sum(a*d for a, d in zip(start, direction, strict=True))
        kinetic = -sum(direction[3:]) * (scale//10)
        assert kinetic >= 0 and all(v <= 0 for v in direction[3:])
        # 先复核较大半径的分支排除，再在共同半径上给出负读数的幅度。
        assert constant + least_sum(weights, bands, branch_radius[axis]) > 0
        margin = constant + least_sum(weights, bands, radius)
        assert margin > arb(".1559")*kinetic
        if kinetic:
            deficits[axis].append((margin/kinetic).lower())
        counts[axis] += 1
    assert list(counts.values()) == [64, 534, 94]
    for axis in BRANCHES:
        print("高层、叶数、分支排除半径、共同半径下亏额下界", axis+1, counts[axis],
              branch_radius[axis].str(5), min(deficits[axis]).str(13), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
