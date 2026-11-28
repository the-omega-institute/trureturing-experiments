"""四枚与七枚已知重谱在二十种谱约定下的共同能带最少段数。"""
from itertools import groupby
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, arb, ctx

TARGETS = (2345712, 2346088, 2346090, 2346094, 2346102, 2346108, 2346122)
LIMITS = ((1, 3), (4, 6), (5, 6), (7, 8), (8, 9), (6, 7), (6, 7))
CHARGES = (1, 2, 3, 4, 6)


def spectra(targets=TARGETS):
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 320
    columns = {}
    for n in targets:
        h = acb.zeta_zero(n).imag
        delta = arb(23)/25 * arb(2)**-20 / (h*h + shell["soft"]["point"]["study"]["offset"])
        phases = shell["geo"]["phases"](h)
        for charge in CHARGES:
            matrix = shell["laplacian"](charge, h, phases=phases)
            eig = matrix.eig(multiple=True)
            assert len(eig) == 19 and all(z.imag.contains(0) and z.real > 0 for z in eig)
            raw = [z.real for z in eig]
            folded = [abs(v * (v - 15 + delta/3)) for v in raw]
            for name, values in (("原谱", raw), ("折叠谱", folded)):
                row = sorted(values, key=lambda v: v.mid())
                assert row[0] > 0 and all(a < b for a, b in zip(row, row[1:]))
                normalized = [v / row[-1] for v in row[:-1]] + [arb(1)]
                for unit, levels in (("原单位", row), ("归一化", normalized)):
                    columns.setdefault((charge, name, unit), []).append(levels)
    return columns


def level_groups(rows):
    events = [(v, i) for i, row in enumerate(rows) for v in row]
    events.sort(key=lambda item: item[0].mid())
    groups = []
    for v, i in events:
        if groups and not groups[-1][0] < v:
            # 各自除以已认证的唯一正最大值后，其比值精确等于 1。
            assert groups[-1][0].is_exact() and v.is_exact() and groups[-1][0] == v == 1
        else:
            groups.append((v, [0] * len(rows)))
        groups[-1][1][i] += 1
    return groups


def optimal_bands(groups, limits):
    # 两个整数数组分别记录最后一个事件组落在能带外、内时的最少段数。
    shape = tuple(hi + 1 for _, hi in limits)
    inf = 200
    assert len(groups) + 1 < inf < 255
    off = np.full(shape, inf, dtype=np.uint8)
    on = np.full_like(off, inf)
    off[(0,) * len(shape)] = 0
    history = []
    for _, inc in groups:
        history.append((off, on))
        src = tuple(slice(0, size - step) for size, step in zip(shape, inc))
        dst = tuple(slice(step, size) for size, step in zip(shape, inc))
        next_on = np.full_like(on, inf)
        next_on[dst] = np.minimum(on[src], off[src] + 1)
        off, on = np.minimum(off, on), next_on

    allowed = tuple(slice(lo, hi + 1) for lo, hi in limits)
    costs = np.minimum(off[allowed], on[allowed])
    optimum = int(costs.min())
    assert optimum < inf
    local = np.unravel_index(costs.argmin(), costs.shape)
    counts = tuple(lo + int(i) for (lo, _), i in zip(limits, local))
    state, cost, inside = counts, optimum, on[counts] == optimum
    selected = []
    for j in reversed(range(len(groups))):
        before_off, before_on = history[j]
        if inside:
            selected.append(j)
            state = tuple(v - d for v, d in zip(state, groups[j][1]))
            assert min(state) >= 0
            inside = before_on[state] == cost
            cost -= int(not inside)
        else:
            inside = before_on[state] == cost
        assert (before_on if inside else before_off)[state] == cost
    assert state == (0,) * len(shape) and cost == 0 and not inside
    selected.reverse()
    runs = []
    for _, consecutive in groupby(enumerate(selected), key=lambda pair: pair[1] - pair[0]):
        indices = [index for _, index in consecutive]
        runs.append((indices[0], indices[-1]))
    assert len(runs) == optimum
    return runs, counts


def check_witness(groups, rows, runs, counts):
    # 每个内部端点取相邻能级外包间隙内的二进有理数，首尾两侧单独处理。
    def gap(j):
        if j == 0:
            return arb(0)
        if j == len(groups):
            return groups[-1][0].upper() + 1
        left, right = groups[j - 1][0], groups[j][0]
        middle = ((left.upper() + right.lower()) / 2).mid()
        assert left < middle < right
        return middle

    intervals = [(gap(a), gap(b + 1)) for a, b in runs]
    for (a, b), (c, d) in zip(intervals, intervals[1:]):
        assert a < b < c < d
    got = tuple(sum(any(a < v <= b for a, b in intervals) for v in row) for row in rows)
    assert got == counts
    return intervals


def check_decimal_example(rows, n):
    # 荷一归一化折叠谱的有理端点，端点单位为百万分之一。
    if n == 4:
        ends = ((915, 43694), (255478, 317055))
        expected = (3, 4, 5, 7)
    else:
        ends = ((135066, 173252), (177694, 252756), (439194, 454774))
        expected = (2, 5, 5, 7, 8, 6, 6)
    bands = [(arb(a)/1000000, arb(b)/1000000) for a, b in ends]
    for row in rows:
        for v in row:
            assert all(v < x or x < v for band in bands for x in band)
    counts = tuple(sum(any(a < v <= b for a, b in bands) for v in row) for row in rows)
    assert counts == expected


def verify():
    three_for_four = {(3, "原谱", "原单位")}
    three_for_seven = three_for_four | {
        (1, "折叠谱", "原单位"), (1, "折叠谱", "归一化"),
        (2, "折叠谱", "原单位"), (4, "原谱", "原单位"),
    }
    for (charge, name, unit), rows in spectra().items():
        for n in (4, 7):
            groups = level_groups(rows[:n])
            assert len(groups) == (19 * n if unit == "原单位" else 18 * n + 1)
            runs, counts = optimal_bands(groups, LIMITS[:n])
            check_witness(groups, rows[:n], runs, counts)
            three_required = three_for_four if n == 4 else three_for_seven
            assert len(runs) == 2 + ((charge, name, unit) in three_required)
            if (charge, name, unit) == (1, "折叠谱", "归一化"):
                check_decimal_example(rows[:n], n)
            ranks = [(a + 1, b + 1) for a, b in runs]
            print(n, charge, name, unit, len(runs), ranks, counts, flush=True)


if __name__ == "__main__":
    verify()
