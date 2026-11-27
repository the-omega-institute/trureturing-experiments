"""仅由实际谱权、门户和单层响应计算的四方向色动能排除判据。"""

from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

DIRECTIONS = ((0, 0, -1, 0), (0, 0, -1, -10), (-20, 29, -1, 0), (-21, 29, -1, -2))
ENERGY = 560000


def projected_tables(core, root, number, portal, *, height=None):
    tables = core["build_grid"](
        root, (0,), number=number, energies=(ENERGY,), portal=portal, height=height
    )
    lower, upper, _, _, start = tables[:5]
    # 端点三式后是 A_QQ/10 和 A_c/10，只保留色坐标。
    coordinates = (0, 1, 2, 4)
    lower = lower[:, coordinates]
    upper = upper[:, coordinates]
    start = [start[i] for i in coordinates]
    ticks = core["cutoff_ticks"]()
    edges = [0, *np.searchsorted(ticks, np.arange(4, 80, 4) * 512), len(ticks) - 1]
    limit = max(int(np.max(np.abs(lower))), int(np.max(np.abs(upper))))
    out = []
    for direction in DIRECTIONS:
        floors, ceilings = core["projection"](lower, upper, direction, limit)
        intervals = [
            (int(floors[a : b + 1].min()), int(ceilings[a : b + 1].max()))
            for a, b in zip(edges[:-1], edges[1:], strict=True)
        ]
        offset = (
            sum(a * d for a, d in zip(start, direction, strict=True)) * core["SCALE"]
        )
        out.append((offset, intervals, core["SCALE"] * sum(map(abs, direction))))
    return out


def support_score(tables, weights, count):
    # 绝对权从大到小；同权时保持原生顺序。已算样本的截断处有严格间隙。
    order = sorted(
        range(len(weights)), key=lambda j: abs(weights[j]).mid(), reverse=True
    )
    major = order[:count]
    rest = [weights[i] for i in order[count:]]
    assert all(w > 0 or w < 0 for w in weights)
    if count:
        assert min(abs(weights[i]).lower() for i in major) > max(
            abs(w).upper() for w in rest
        )
    p = sum((w for w in rest if w > 0), arb(0))
    n = -sum((w for w in rest if w < 0), arb(0))
    projected = []
    for offset, spans, scale in tables:
        neutral = offset + p * min(v[0] for v in spans) - n * max(v[1] for v in spans)
        choices = [
            [weights[i] * (l if weights[i] > 0 else h) / scale for l, h in spans]
            for i in major
        ]
        projected.append((neutral / scale, choices))
    lowest = None
    for cell in product(range(len(tables[0][1])), repeat=count):
        values = [
            base + sum((d[k] for d, k in zip(choices, cell, strict=True)), arb(0))
            for base, choices in projected
        ]
        best = values[0]
        for value in values[1:]:
            best = best.max(value)
        lowest = best if lowest is None else lowest.min(best)
    return lowest


def evaluate(root, number, *, portal=None):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    tables = projected_tables(core, root, number, portal)
    spectra = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    rows = next(v for v in spectra["shell_weights"](root, number) if v[0] == 6)[2]
    result = []
    for row in rows:
        weights = [w.real for w in row]
        assert sum(weights).contains(1)
        result.append([support_score(tables, weights, k) for k in (0, 2, 3)])
    return result


def verify(root):
    for shell in range(4, 8):
        scores = evaluate(root, 2_306_888_100_085 + shell)
        for j, (none, two, three) in enumerate(scores, 1):
            print(
                "壳、枢轴、分别保留零/二/三个谱权的指标",
                shell,
                j,
                *[v.str(13) for v in (none, two, three)],
                flush=True,
            )
        assert scores[0][0] > arb(".17")
        assert scores[0][2] > arb(".44")
        assert scores[1][2] > arb(".10") and scores[1][1] < arb("-2.98")
        assert scores[2][2] < arb("-24")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
