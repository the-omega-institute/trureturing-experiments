"""由壳上给定的坐标三循环区分固定见证的四个二次源协方差。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify():
    neighbours = run_path(str(Path(__file__).with_name("ch39-four-quadratic-source-neighbours.py")))
    load = neighbours["load"]
    model, _, _ = load("two-zero-orbit-condition")["calibration"]()
    gram = load("quadratic-gram-source")
    adj = model["adj"]
    sites = model["shell"]["geo"]["sites"]
    permutation = [sites.index((y, z, x)) for x, y, z in sites]
    rotation = acb_mat([[int(i == j) for j in permutation] for i in range(len(sites))])
    source = model["sources"]
    observable = acb(0, 1) * adj(source) * (rotation - adj(rotation)) * source
    a, b = [gram["at_zero"](model, q, index) for q, index in ((-4, 2345666), (2, 964890))]
    words = neighbours["words"](adj(a)*a, adj(b)*b, gram)
    chosen = list(neighbours["coefficients_at_roots"](model, a, b))
    results = []
    for (center, coeff), expected in zip(chosen, ("3.077473", "3.273211", "3.240416", "3.161874"), strict=True):
        covariance = sum((v*w for v, w in zip(coeff, words, strict=True)), acb_mat(3, 3))
        value = (observable*covariance).trace().real
        assert abs(value - arb(expected)) < arb(".00001")
        results.append(value)
        print("第一相位盒中心与旋转读数", center[0].str(7), value.str(9))
    assert all(abs(a-b) > arb(".0327") for a, b in combinations(results, 2))


if __name__ == "__main__":
    verify()
