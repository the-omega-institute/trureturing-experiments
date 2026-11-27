"""第八壳三条贝伦德最短弧的相干站位读出与四方向支撑判据。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, arb


SCALE = 10**7
NUMBER = 2_306_888_100_093


def polynomials(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    near = shell["boundary_sites"]()
    pairs = [(sum((a-b)**2 for a, b in zip(sites[i], sites[j], strict=True)), i, j)
             for i, j in combinations(range(19), 2)]
    bridges = [(i, j) for chord, i, j in pairs if chord == min(t[0] for t in pairs)]
    assert bridges == [(6, 11), (7, 12), (8, 18)] and near == [8, 7, 6]
    far = [dict(bridges)[i] for i in near]
    assert all(sum(a*a for a in site) == 867 for site in sites)
    for i, j in bridges:
        midpoint = Fraction(sum((a+b)**2 for a, b in zip(sites[i], sites[j], strict=True)), 4)
        assert midpoint == 864 and midpoint / 867 == Fraction(288, 289)
    _, vectors = lib("adjacent-pvm-transport")["spectrum"](shell, 6, NUMBER)
    matrix = shell["laplacian"](6, acb.zeta_zero(NUMBER).imag)
    links = [-matrix[a, b] for a, b in zip(near, far, strict=True)]
    assert all(z.contains(1) for z in links)
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    pole = shell["soft"]["pole"]
    output = []
    single_edge = [[] for _ in near]
    t = 3 * arb.pi() / 8
    for pivot in range(1, 4):
        rows = []
        weights = gram[pivot]
        for j in range(19):
            left = [vectors[a, j] for a in near]
            right = [link * vectors[b, j] for link, b in zip(links, far, strict=True)]
            if pivot == 2:
                for which in range(3):
                    measured = [*left]
                    measured[which] = left[which]*t.cos() + right[which]*t.sin()
                    source = weights[0]*vectors[pole, j] + sum(
                        c*z for c, z in zip(weights[1:], measured, strict=True)
                    )
                    single_edge[which].append((source*measured[pivot-1].conjugate()/weights[pivot]).real)
            aa = sum(c*z for c, z in zip(weights[1:], left, strict=True))
            ab = sum(c*z for c, z in zip(weights[1:], right, strict=True))
            p = weights[0] * vectors[pole, j]
            q = left[pivot-1].conjugate() / weights[pivot]
            r = right[pivot-1].conjugate() / weights[pivot]
            # 由 cos²、sin²、cos·sin 换成双角，避免在同一参数区间重复取宽。
            rows.append(tuple(z.real for z in (
                (aa*q+ab*r)/2, p*q, p*r, (aa*q-ab*r)/2, (aa*r+ab*q)/2
            )))
        output.append(rows)
    assert all(sum(row).contains(1) for row in single_edge)
    return output, single_edge


def angle_weights(rows, s):
    t = arb.pi() * s / 2
    harmonics = (1, t.cos(), t.sin(), (2*t).cos(), (2*t).sin())
    return [sum(x*y for x, y in zip(row, harmonics, strict=True)) for row in rows]


def integer_bounds(value):
    scaled = value * SCALE
    return int(scaled.lower().floor().unique_fmpz()), int(scaled.upper().ceil().unique_fmpz())


def score_box(tables, weights, major):
    rest = [w for i, w in enumerate(weights) if i not in major]
    pos = sum(w.max(arb(0)) for w in rest)
    neg = sum(w.min(arb(0)) for w in rest)
    projections = []
    for offset, bands, scale in tables:
        base = integer_bounds((offset + pos*min(t[0] for t in bands)
                               + neg*max(t[1] for t in bands)) / scale)
        major_terms = [
            [integer_bounds((weights[i]*l).min(weights[i]*h) / scale) for l, h in bands]
            for i in major
        ]
        # 四个整数相加，下界、上界分别在最后一维；预先排除 int64 溢出。
        assert max(abs(z) for pair in [base, *(t for col in major_terms for t in col)]
                   for z in pair) < 2**59
        a, b, c = np.asarray(major_terms, dtype=np.int64)
        projections.append(base + a[:, None, None, :] + b[None, :, None, :] + c[None, None, :, :])
    return np.max(projections, axis=0).min(axis=(0, 1, 2))


def score(tables, weights):
    critical = sorted((abs(w).lower() for w in weights), reverse=True)[2]
    candidates = [i for i, w in enumerate(weights) if not abs(w).upper() < critical]
    if len(candidates) > 5:
        return None
    bounds = []
    for major in combinations(candidates, 3):
        # 三层之外若有必然更大的权，则该组合不可能入选。
        outside = max(abs(w).lower() for i, w in enumerate(weights) if i not in major)
        if all(not abs(weights[i]).upper() < outside for i in major):
            bounds.append(score_box(tables, weights, major))
    assert bounds
    return min(a for a, _ in bounds), max(b for _, b in bounds)


def cover(tables, rows, start, stop, greater, limit):
    stack = [(arb(start), arb(stop), 0)]
    accepted = []
    while stack:
        l, r, depth = stack.pop()
        bounds = score(tables, angle_weights(rows, l.union(r)))
        if bounds is not None:
            lo, hi = bounds
            if (greater and lo > limit) or (not greater and hi < limit):
                accepted.append(lo if greater else hi)
                continue
        assert depth < 22, (l, r, bounds)
        mid = (l+r)/2
        stack.extend(((l, mid, depth+1), (mid, r, depth+1)))
    verified = min(accepted) if greater else max(accepted)
    print("弧分数、闭段数、指标有向界", start, stop, len(accepted), verified / SCALE, flush=True)


def verify(root):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    support = run_path(str(root / "ch39-cross-shell-pivot-support.py"))
    tables = support["projected_tables"](core, root, NUMBER, None)
    rows, single_edge = polynomials(root)
    for which, (weights, predict_positive) in enumerate(zip(single_edge, (True, False, True), strict=True), 1):
        lo, hi = score(tables, weights)
        outcome = ("命中" if predict_positive else "失配") if lo > 0 else (
            ("失配" if predict_positive else "命中") if hi < 0 else "未定"
        )
        print("只动一条最短边、弧分数 3/4、第二枢轴支撑包围、预登记结果",
              which, int(lo)/SCALE, int(hi)/SCALE, outcome, flush=True)
    for s in (0, 1):
        result = [score(tables, angle_weights(row, arb(s))) for row in rows]
        print("站位端、三个指标的整数包围（单位一千万分之一）", s, result, flush=True)
    for pivot, left, right, greater, limit in (
        (1, "0", "1", True, 3500000),
        (2, "0", ".62", True, 500000),
        (2, ".71", ".80", False, -2000000),
        (2, ".84", "1", True, 1500000),
        (3, "0", "1", False, -10000000),
    ):
        print("枢轴", pivot, flush=True)
        cover(tables, rows[pivot-1], left, right, greater, limit)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
