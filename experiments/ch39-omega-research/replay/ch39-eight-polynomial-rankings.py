"""用整数恒等式复核 32 个七壳最小支撑的严格多项式排序次数。"""
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
from runpy import run_path


def verify(root, winners):
    points = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))["geo"]["sites"]
    powers = [term for degree in range(4)
              for term in combinations_with_replacement(range(3), degree)]
    design = [[prod(point[axis]-17 for axis in term) for term in powers] for point in points]
    records = json.loads((root / "ch39-eight-polynomial-rankings.json").read_text())
    assert len(records) == 32 and {tuple(s) for s, _, _ in records} == set(winners)
    quadratic = 0
    for support, coeff, dual in records:
        size = 10 if dual is None else 20
        assert len(coeff) == size and all(type(c) is int for c in coeff)
        signed = [[(1 if i in support else -1)*x for x in row] for i, row in enumerate(design)]
        assert all(sum(a*b for a, b in zip(row[:size], coeff, strict=True)) > 0 for row in signed)
        if dual is None:
            quadratic += 1
            continue
        # 正整数组合湮灭全部二次单项式（含常数），故不可能存在严格分隔。
        assert dual and all(0 <= i < 19 and type(w) is int and w > 0 for i, w in dual)
        assert all(sum(w*signed[i][j] for i,w in dual) == 0 for j in range(10))
    assert quadratic == 11
    for axis, support in ((1,(1,2,3,4,7,10,12,15)), (2,(0,2,3,5,6,9,11,13))):
        scores = [(site[axis]-13)**2 for site in points]
        assert support in winners
        assert min(scores[i] for i in support) == 100
        assert max(scores[i] for i in range(19) if i not in support) == 64
    print("严格排序所需最低次数为二、三的支撑数", quadratic, len(winners)-quadratic)
    print("两个坐标平方评分的内外边界", 100, 64)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    classify = run_path(str(root / "ch39-eight-optimal-classification.py"))["experiment"]
    verify(root, classify(root, *source["verify"](root, return_operator="spectral")))
