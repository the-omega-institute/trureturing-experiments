"""不含原壳副本且至多 217 点时，单钟、共参考双钟、共站双钟的最小误差。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def nearest_tip(angles, target, index):
    # 基准角严格递增，两个目标角的插入位置各自只会前移。
    while index < len(angles) and angles[index] < target:
        index += 1
    if index < len(angles):
        assert angles[index] > target
    neighbors = range(max(index-1, 0), min(index+1, len(angles)))
    errors = [(abs(angles[j]-target), j) for j in neighbors]
    least, at = min(errors, key=lambda entry: float(entry[0]))
    assert all(j == at or least < e for e, j in errors)
    return index, at, least


def test_list(lengths, u, frequencies, records, site=None):
    angles = [(1-arb(d)/(6*u*u)).acos() for d in lengths]
    assert all(a < b for a, b in zip(angles, angles[1:]))
    indices = [0, 0]
    for ref, base in enumerate(angles):
        tips, errors = [], []
        for k, omega in enumerate(frequencies):
            indices[k], tip, err = nearest_tip(angles, base+omega, indices[k])
            tips.append(lengths[tip])
            errors.append(err)
        if site is None:
            for k, err in enumerate(errors):
                records[k].append((err, u, (lengths[ref], tips[k]), None))
        key = 2 if site is None else 3
        records[key].append((errors[0].max(errors[1]), u, (lengths[ref], *tips), site))


def verify():
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-nineteen-lattice-shells.py"))
    ctx.prec = 144
    frequencies = [(arb(b)/867).acos()-(arb(745)/867).acos() for b in (505, 465)]
    records = [[], [], [], []]
    scanned = 0
    distance = source["square_distance"]
    for u in range(1, 218, 2):
        if u % 17 == 0:
            continue
        points = source["shell"](3*u*u)
        if len(points) > 217:
            continue
        scanned += 1
        lengths = sorted({distance(p, q) for p, q in combinations(points, 2)} | {0})
        test_list(lengths, u, frequencies, records)
        for site in sorted({tuple(sorted(p)) for p in points}):
            test_list(sorted({distance(site, q) for q in points}), u, frequencies, records, site)
    assert scanned == 87
    expected = [(211, (41428,117556), None), (157, (216,19684), None),
                (211, (41428,117556,130104), None), (105, (196,7756,9576), (35,119,133))]
    for rows, answer in zip(records, expected, strict=True):
        row = min(rows, key=lambda entry: float(entry[0]))
        best, *key = row
        assert tuple(key) == answer
        assert all(candidate is row or best < candidate[0] for candidate in rows)
        print("最优角频误差、奇因子、参考及钟端的平方弦长、公共站", best.str(16), *key)


if __name__ == "__main__":
    verify()
