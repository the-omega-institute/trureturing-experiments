"""用八维实 Gram 表示逐项重建十九点壳在截止 216 的三时钟矩条件。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import fmpq
from flint.types._gr import gr_real_qqbar_ctx

plane = run_path(str(Path(__file__).with_name("ch39-three-clock-plane-completion.py")))


def verify():
    real = gr_real_qqbar_ctx.new()
    zero = real(0)
    points, distance, g0, delta, entries = plane["prepare"](fmpq(1))

    def outer(i, j):
        a, b, _, _ = entries(i, j)
        return real(-a / b)

    def line(i, j, z):
        a, b, c, d = map(real, entries(i, j))
        return a + b*z, c + d*z

    def residual(i, j, left, right):
        a, b = line(i, j, left)
        return a + b*right

    x, w = outer(4, 2), outer(8, 0)
    A, B = [residual(4, j, x, x) for j in (12, 15)]
    v_squared = (2*A - B)/3
    assert real.lt(zero, v_squared) and real.lt(zero, 2*A + B)
    # 平面坐标用 v 作单位，q=u/v，内积的权重为 v²。
    q = (((2*A + B)/v_squared).sqrt() - 1)/2
    a7, b7 = line(4, 7, x)
    a16, b16 = line(4, 16, x)
    y = (a16 - q*a7)/(q*b7 - b16)
    r = (a7 + b7*y)/(3*v_squared)
    s = residual(4, 14, x, w)/(3*v_squared*q)

    # 这条正交直线在点 6 的轨道取单位坐标，在点 8 的轨道取 k。
    weight = 3*v_squared*r*r - residual(6, 7, y, y)
    assert real.lt(zero, weight)
    k = (3*v_squared*r*s - residual(6, 9, y, w))/weight
    rows = {}
    for i, z, vector, scalar in [(4, x, (q, real(1), -q-1), zero),
                                 (6, y, (r, r, -2*r), real(1)),
                                 (8, w, (s, s, -2*s), k)]:
        for perm in permutations(range(3)):
            j = points.index(tuple(points[i][h] for h in perm))
            row = (z, tuple(vector[h] for h in perm), scalar)
            if j in rows:
                assert rows[j] == row
            rows[j] = row

    for i in (*plane["outer"], plane["centre"]):
        rows[i] = (zero, (zero, zero, zero), zero)
    assert len(rows) == 19
    trace = zero
    fixed = 0
    for i in range(19):
        zi, vi, si = rows[i]
        for j in range(i, 19):
            if i != j and distance[i][j] <= 216:
                continue
            zj, vj, sj = rows[j]
            diff = residual(i, j, zi, zj)
            diff += v_squared*sum(a*b for a, b in zip(vi, vj)) + weight*si*sj
            if i == j:
                trace += (real(delta[i, i]) + diff)/10**10
            else:
                assert real.is_zero(diff), (i, j)
                fixed += 1
    inv_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    assert fixed == 114
    assert real.lt(trace*real(inv_norm), real((1-fmpq(98, 1000))**2))
    assert real.lt(trace, real(fmpq(2, 10**6)))
    print("Cutoff 216, at most 27 auxiliary states,", fixed,
          "exact prescribed edges, Gram trace < 2 * 10^-6")


if __name__ == "__main__":
    verify()
