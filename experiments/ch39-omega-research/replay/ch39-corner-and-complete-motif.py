"""穷举同端点弦比、三边相似图形与完整十九点壳的嵌入。"""

from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sys


def find_corners(points, distance):
    edges = defaultdict(list)
    for p, q in combinations(points, 2):
        edges[distance(p, q)].append((p, q))
    corners = []
    complete = []
    for b, short_edges in edges.items():
        if b % 8 or 181 * (b // 8) not in edges:
            continue
        a = 181 * (b // 8)
        neighbors = defaultdict(list)
        for p, q in edges[a]:
            neighbors[p].append(q)
            neighbors[q].append(p)
        for p, q in short_edges:
            for apex, short_end in ((p, q), (q, p)):
                for long_end in neighbors[apex]:
                    c = distance(short_end, long_end)
                    item = (apex, short_end, long_end, b, a, c)
                    corners.append(item)
                    if 8 * c == 201 * b:
                        complete.append(item)
    return corners, complete


def find_copies(u, points, reference, distance):
    m = u // 17
    anchors = ((17, 17, 17), (13, 13, 23), (13, 23, 13))
    # 三个锚向量线性无关。同一原点球上，到它们的平方距离唯一决定一点。
    keys = {tuple(m**2 * distance(p, a) for a in anchors) for p in reference}
    assert len(keys) == len(reference) == 19
    table = {p: {q: distance(p, q) for q in points} for p in points}
    near_length = 68 * m**2
    other_length = 200 * m**2
    candidates = 0
    copies = Counter()
    for a in points:
        near = [p for p in points if table[a][p] == near_length]
        for b, c in permutations(near, 2):
            if table[b][c] != other_length:
                continue
            candidates += 1
            image = frozenset(p for p in points if (table[a][p], table[b][p], table[c][p]) in keys)
            if len(image) == len(reference):
                copies[image] += 1
    standard = frozenset(tuple(m * x for x in p) for p in reference)
    assert copies == {standard: 6}
    return candidates, standard


def verify(root):
    lattice = run_path(str(root / "ch39-nineteen-lattice-shells.py"))
    shell, distance = lattice["shell"], lattice["square_distance"]
    reference = shell(867)
    assert abs(lattice["volume6"]((0, 0, 0), ((17, 17, 17), (13, 13, 23), (13, 23, 13)))) == 1700
    # N(u)≥u 说明此奇数列表覆盖所有不超过 217 点的完整对角壳相似类。
    bound = 217
    checked = Counter()
    ties = defaultdict(list)
    corner_minima = []
    complete_minima = []
    copies_report = []
    for u in range(1, bound + 1, 2):
        points = shell(3 * u**2)
        n = len(points)
        assert n >= u
        if n > bound:
            continue
        checked[u % 17 != 0] += 1
        corners, complete = find_corners(points, distance)
        if n in (49, 217):
            ties[n].append((u, len(corners), len(complete)))
        if u % 17:
            if n <= 49 and corners:
                corner_minima.append((n, u, len(corners)))
                assert set(item[3:] for item in corners) == {(416, 9412, 8436)}
                print("最小共端点弦例（顶点、短端、长端、三边平方）", min(corners))
            if complete:
                complete_minima.append((n, u, len(complete)))
                assert set(item[3:] for item in complete) == {(2336, 52852, 58692)}
                print("最小完整三角例（顶点、短端、长端、三边平方）", min(complete))
        else:
            candidates, standard = find_copies(u, points, reference, distance)
            outside = [item for item in complete if not set(item[:3]) <= standard]
            copies_report.append((u, n, candidates, len(corners), len(complete), len(outside)))
            if outside:
                print("唯一整壳副本外的完整相似三角", u, min(outside))
    assert checked == {True: 87, False: 4}
    assert corner_minima == [(49, 49, 6)]
    assert complete_minima == [(217, 217, 12)]
    assert ties == {49: [(35, 0, 0), (47, 0, 0), (49, 6, 0)],
                    217: [(145, 18, 0), (155, 18, 0), (203, 0, 0), (217, 18, 12)]}
    assert copies_report == [(17, 19, 6, 6, 6, 0), (51, 76, 12, 12, 6, 0),
                             (85, 133, 18, 36, 12, 6), (119, 133, 18, 30, 12, 6)]
    print("候选壳数（17 不整除、17 整除）", checked[True], checked[False])
    print("两个临界点数各壳的奇因子、共端点弦对数、完整三角数", dict(ties))
    print("含副本壳的奇因子、点数、候选锚像、弦对数、完整三角数、整副本外的完整三角数",
          copies_report)


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)
