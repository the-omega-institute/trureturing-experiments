"""保持球面和坐标置换对称的形变下，无量纲几何标尺的变化。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path


def verify():
    root = Path(__file__).resolve().parent
    geometry = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))["graph"]
    sites = geometry["sites"]
    assert all(sum(v*v for v in site) == 867 for site in sites)
    heights = [Fraction(sum(site), 51) for site in sites]
    chords = {(i, j): sum((a-b)**2 for a, b in zip(sites[i], sites[j], strict=True))
              for i, j in combinations(range(len(sites)), 2)}
    cap = Fraction(1, 500)
    # 不同原始弦距的次序在全形变闭段上均严格保持。
    levels = sorted(set(chords.values()))
    min_ratio = min(Fraction(hi, lo) for lo, hi in zip(levels, levels[1:]))
    assert min_ratio > ((1 + cap) / (1 - cap))**2
    for distance, sums in ((32, (35, 35)), (724, (35, 49)), (804, (35, 49))):
        assert {tuple(sorted((sum(sites[i]), sum(sites[j]))))
                for (i, j), d in chords.items() if d == distance} == {sums}
    # 距中心轴上固定点的最近、最远层，以及三种决定标尺的弦距，均不分裂。
    pole = sites.index((17, 17, 17))
    by_radius = {}
    for i, site in enumerate(sites):
        if i != pole:
            radius = chords[tuple(sorted((pole, i)))]
            by_radius.setdefault(radius, set()).add(heights[i])
    assert all(len(group) == 1 for group in by_radius.values())
    def scale(s):
        return 724 * (51 + 35*s) / (51 + 49*s)
    assert 723 < scale(cap) < 724 < scale(-cap) < 725
    print("最小相邻平方弦距比与最大相对次序扰动", min_ratio, ((1+cap)/(1-cap))**2)
    print("两端标尺、中心导数", scale(cap), scale(-cap), Fraction(-724*14, 51))


if __name__ == "__main__":
    verify()
