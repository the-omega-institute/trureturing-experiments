"""核对两只相位钟的共同基点、独立平方类和有限名单的闭路长度。"""
from fractions import Fraction
from itertools import combinations
from math import isqrt, prod
from pathlib import Path
from runpy import run_path

root = Path(__file__).parent


def verify():
    geometry = run_path(str(root / "ch39-arc-relative-determinant.py"))
    short = run_path(str(root / "ch39-short-arc-born-weights.py"))
    pole = geometry["sites"].index((17, 17, 17))
    lengths = {(i, j): length for i, j, length, _ in geometry["edges"]}
    for mobile in (724, 804):
        i, j = next(edge for edge, length in lengths.items() if length == mobile)
        spokes = [lengths[tuple(sorted((pole, v)))] for v in (i, j)]
        assert all(length not in (724, 804) for length in spokes)
        print("极点闭路、活动边、两条辐边", (pole, i, j, pole), mobile, spokes)

    radius_squared = geometry["R"]
    lengths_used = (244, 724, 804)
    radicands = [-d * (4 * radius_squared - d) for d in lengths_used]
    for size in (1, 2, 3):
        for subset in combinations(radicands, size):
            value = prod(subset)
            assert value < 0 or isqrt(value)**2 != value
    double_cosine = 2 - Fraction(724, radius_squared)
    assert double_cosine.denominator > 1
    print("三组独立平方类，活动角的二倍余弦", double_cosine)

    # 对应已归档的八列正解，整个十五点约束名单也加入插值。
    targets = run_path(str(root / "ch39-five-neutral-quark-probability.py"))["eight_indices"]
    count = len(set(short["indices"]) | set(targets))
    assert count == 22
    print("有限名单编号数、插值使用的最大边数", count, 3*(count-1))


if __name__ == "__main__":
    verify()
