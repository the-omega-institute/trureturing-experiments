"""直接取有理圆弧读出，定位第八壳第二枢轴的两个变号及其前三层。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb


def interval(value):
    q = Fraction(value)
    return arb(q.numerator) / q.denominator


def rational_weights(rows, t, *, derivative=False):
    d = 1+t*t
    u, v = (1-t*t)/d, 2*t/d
    terms = (1, u, v, u*u-v*v, 2*u*v)
    if derivative:
        du, dv = -4*t/d**2, 2*(1-t*t)/d**2
        terms = (0, du, dv, 2*(u*du-v*dv), 2*(du*v+u*dv))
    return [sum(a*b for a, b in zip(row, terms, strict=True)) for row in rows]


def check_sign(arc, tables, rows, left, right, sign):
    stack = [(interval(left), interval(right), 0)]
    count = 0
    while stack:
        a, b, depth = stack.pop()
        bounds = arc["score"](tables, rational_weights(rows, a.union(b)))
        if bounds is not None and (bounds[0] > 0 if sign > 0 else bounds[1] < 0):
            count += 1
            continue
        assert depth < 26, (a, b, bounds)
        mid = (a+b)/2
        stack.extend(((a, mid, depth+1), (mid, b, depth+1)))
    print("有理区间、符号、闭段数", left, right, sign, count, flush=True)


def crossing(arc, tables, rows, l, r, major, direction):
    box = interval(l).union(interval(r))
    weights = rational_weights(rows, box)
    assert all(w > 0 or w < 0 for w in weights)
    assert min(abs(weights[i]).lower() for i in major) > max(
        abs(w).upper() for i, w in enumerate(weights) if i not in major
    )
    zero_offsets = [(arb(0), bands, scale) for _, bands, scale in tables]
    value = arc["projected_box"](tables, weights, major)
    slope = arc["projected_box"](
        zero_offsets, rational_weights(rows, box, derivative=True), major,
        signs=[w > 0 for w in weights],
    )
    # 可能取得外层最小值的截止盒及其中可能取到最大值的方向。
    lower, upper = value[..., 0], value[..., 1]
    cell = lower.max(axis=0) <= upper.max(axis=0).min()
    possible = cell[None, ...] & (upper >= lower.max(axis=0)[None, ...])
    assert np.any(possible)
    bound = (slope[..., 0][possible].min() if direction > 0 else slope[..., 1][possible].max())
    assert direction * bound > 0
    print("变号窗、固定前三层、各可能活动分支的有向导数界",
          l, r, tuple(i+1 for i in major), int(bound)/arc["SCALE"], flush=True)


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    tables = lib("cross-shell-pivot-support")["projected_tables"](core, root, arc["NUMBER"], None)
    all_rows, _ = arc["polynomials"](root, include_single_edge=False)
    rows = all_rows[1]
    for a, b, c, sign in ((3,4,5,1), (5,12,13,-1), (7,24,25,1)):
        assert a*a+b*b == c*c
        t = Fraction(b, c+a)
        weights = rational_weights(rows, interval(t))
        assert sum(weights).contains(1)
        lo, hi = arc["score"](tables, weights)
        assert (lo > 0 if sign > 0 else hi < 0)
        print("勾股单位向量、半角参数、第二枢轴包围", (a,b,c), t,
              int(lo)/arc["SCALE"], int(hi)/arc["SCALE"], flush=True)
    for l, r, sign in (("0","0.60328",1), ("0.60330","0.73615",-1), ("0.73617","1",1)):
        check_sign(arc, tables, rows, l, r, sign)
    crossing(arc, tables, rows, "0.60328", "0.60330", (0,15,16), -1)
    crossing(arc, tables, rows, "0.73615", "0.73617", (0,7,16), 1)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
