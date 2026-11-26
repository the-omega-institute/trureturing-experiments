"""完整对角整数壳中，十九点以后首次出现平方弦长比 181/8 的点数。"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys


def verify(root):
    lattice = run_path(str(root / "ch39-nineteen-lattice-shells.py"))
    shell = lattice["shell"]
    distance = lattice["square_distance"]
    # 已有 37 点见证。类数局部因子逐个至少为 p^k，故点数至多 37 时 u 至多 37。
    upper = 37
    scanned = defaultdict(list)
    hits = []
    representatives = {
        17: (((1, 5, 29), (13, 23, 13)), ((1, 5, 29), (5, 1, 29))),
        25: (((5, 13, 41), (17, 31, 25)), ((1, 5, 43), (5, 1, 43))),
        37: (((11, 31, 55), (23, 47, 37)), ((7, 43, 47), (7, 47, 43))),
    }
    for u in range(1, upper + 1, 2):
        points = shell(3 * u**2)
        count = len(points)
        assert count >= u
        if count > upper:
            continue
        scanned[count].append(u)
        chords = defaultdict(list)
        for edge in combinations(points, 2):
            chords[distance(*edge)].append(edge)
        length_pairs = sorted((181 * (b // 8), b) for b in chords
                              if b % 8 == 0 and 181 * (b // 8) in chords)
        if not length_pairs:
            continue
        hits.append((count, u))
        assert length_pairs == [(724, 32)]
        large, small = chords[724], chords[32]
        assert tuple(map(len, (large, small))) == (6, 3)
        assert all(e in options for e, options in zip(representatives[u], (large, small)))
        third_sides = Counter()
        for e, f in product(large, small):
            endpoints = set(e) ^ set(f)
            if len(endpoints) == 2:
                third_sides[distance(*endpoints)] += 1
        assert third_sides == ({804: 6} if u == 17 else {})
        assert len(chords) == {17: 30, 25: 89, 37: 92}[u]
        print("点数、奇因子、半径平方、不同平方弦长数", count, u, 3 * u**2, len(chords))
        print("724 弦与 32 弦各一例", representatives[u])
        print("两类弦数、共端点时第三边平方及重数", (len(large), len(small)), dict(third_sides))
    assert dict(scanned) == {
        1: [1], 4: [3], 7: [5, 7], 13: [9, 11, 13], 19: [17, 19], 25: [23],
        28: [15, 21], 31: [29, 31], 37: [25, 37],
    }
    assert hits == [(19, 17), (37, 25), (37, 37)]
    print("点数至多 37 的全部奇因子，按点数分组", dict(sorted(scanned.items())))
    # 任意完整相似副本与母壳同外接球。其 32 弦平方变成 32 (u/17)^2。
    assert all(Fraction(32 * u**2, 17**2).denominator == 17**2 for u in (25, 37))
    original = shell(867)
    larger = shell(3 * 51**2)
    assert len(larger) == 76 and {tuple(3 * x for x in p) for p in original} <= set(larger)
    print("包含 3 倍原壳的完整壳，奇因子与总点数", 51, len(larger))


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)
