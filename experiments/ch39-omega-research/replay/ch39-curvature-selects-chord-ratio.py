"""曲率必要窗在四种轨道数与全部原点对弦距比中仅保留一组。"""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from runpy import run_path


def verify():
    root = Path(__file__).resolve().parent
    sites = run_path(str(root / "ch39-arc-relative-determinant.py"))["sites"]
    geometry = run_path(str(root / "ch39-scale-free-geometric-pair.py"))
    _, _, _, fixed = geometry["readout"](sites)
    active, = [count for _, _, count in fixed if count]
    assert active == 96 and len(fixed) == 6
    chords = sorted({sum((a-b)**2 for a, b in zip(x, y, strict=True))
                     for x, y in combinations(sites, 2)})
    assert len(chords) == 30
    lower, upper = Fraction("31.9830"), Fraction("32.0170")
    passing = []
    below, above = [], []
    for order, a, b in product((1, 2, 3, 6), chords, chords):
        k = active // order
        d = k * Fraction(a, b)
        effective = k * (724 / d)**2
        record = effective, k, a, b
        if effective < lower:
            below.append(record)
        elif effective > upper:
            above.append(record)
        else:
            passing.append(record)
    assert passing == [(32, 32, 724, 32)]
    assert len(below) + len(above) + len(passing) == 3600
    assert max(below)[0] == Fraction(98283, 3136)
    assert min(above)[0] == Fraction(32761, 1014)
    print("合格的等效倍率、轨道数、两个平方弦距", passing)
    print("两侧最近的排除值", max(below)[0], min(above)[0])


if __name__ == "__main__":
    verify()
