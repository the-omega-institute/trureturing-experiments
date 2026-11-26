"""若 Gram 补全保持三坐标循环置换，截止低于 392 时排除五维。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import fmpq_poly

pairs = run_path(str(Path(__file__).with_name("ch39-three-clock-outer-pairs.py")))
previous = pairs["previous"]


def determinant(m):
    if len(m) == 1:
        return m[0][0]
    return sum((-1)**k * x * determinant([row[:k] + row[k+1:] for row in m[1:]])
               for k, x in enumerate(m[0]))


def positive_for_positive_t(p):
    assert p and all(c >= 0 for c in p.coeffs())


def verify():
    _, _, _, g0, g1, g2 = previous["MEM"]["coefficients"]()
    delta = g2 - g1 * g0.solve(g1)
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17, 17, 17)})
    outer = pairs["outer"]
    t = fmpq_poly([0, 1])
    g0_outer = pairs["projected_outer_gram"](delta, 0)
    g1_outer = pairs["projected_outer_gram"](delta, 1)
    g = [[g0_outer[i, j] + t * (g1_outer[i, j] - g0_outer[i, j])
          for j in range(6)] for i in range(6)]
    # K(t) 对 t 为仿射。K(0) 正半定，斜率在成对常向量上是三点图拉普拉斯的倍数。
    for i in range(6):
        for j in range(6):
            assert g[i][j][1] == (2 if i // 2 == j // 2 else -1)
    signs = (1, -1, -1, 1, 1, -1)
    assert all(sum(row) == 0 and sum(s * p for s, p in zip(signs, row)) == 0 for row in g)
    orbit = set()
    for perm in permutations(range(3)):
        idx = [points.index(tuple(p[k] for k in perm)) for p in points]
        lookup = [outer.index(idx[i]) for i in outer]
        assert all(g[i][j] == g[lookup[i]][lookup[j]] for i in range(6) for j in range(6))
        assert all(delta[i, j] == delta[idx[i], idx[j]] for i in range(19) for j in range(19))
        orbit.add(idx[4])
    assert {4, 12, 15} <= orbit
    def distance(i, j):
        return local["squared"](points[i], points[j])

    def lifted(i):
        fixed = [k for k, j in enumerate(outer) if distance(i, j) >= 392]
        assert len(fixed) == 4
        m = [[g[k][l] for l in fixed] for k in fixed]
        den = determinant(m)
        positive_for_positive_t(den)
        coeff = [[fmpq_poly() for _ in range(7)] for _ in range(2)]
        for j, k in enumerate(fixed):
            for degree in range(2):
                # 克拉默法则。参数 x 是该点在中心方向上的无量纲坐标。
                rhs = [-t if degree else delta[i, outer[l]] for l in fixed]
                replaced = [row[:j] + [rhs[l]] + row[j+1:] for l, row in enumerate(m)]
                coeff[degree][k] = determinant(replaced)
        coeff[1][6] = den
        return den, coeff

    den, base = lifted(4)
    def equation(j):
        other_den, other = lifted(j)
        assert other_den == den and distance(4, j) >= 392
        def inner(x, y):
            return (sum(x[k] * g[k][l] * y[l] for k in range(6) for l in range(6))
                    + t * x[6] * y[6])
        h = [[inner(x, y) for y in other] for x in base]
        h[0][0] -= delta[4, j] * den**2
        return h

    h, k = equation(12), equation(15)
    anti = h[0][1] - h[1][0]
    def constant_sign(p):
        positive_for_positive_t(p if p.leading_coefficient() > 0 else -p)
    constant_sign(anti)
    # 三条有向的内点边依次给同一反分式线性变换，内点坐标可各自任意。
    det_h = h[0][0]*h[1][1] - h[0][1]*h[1][0]
    positive_for_positive_t(det_h - anti**2)
    rotate = [points.index((p[1],p[2],p[0])) for p in points]
    cyclic = (4,12,13)
    assert ([rotate[i] for i in cyclic] in ([12,13,4],[13,4,12])
            and all(distance(cyclic[i],cyclic[(i+1)%3])==392 for i in range(3)))
    assert h[1][1] == k[1][1]
    slope = h[0][1] + h[1][0] - k[0][1] - k[1][0]
    intercept = h[0][0] - k[0][0]
    constant_sign(slope)
    obstruction = k[0][0]*slope**2 - (k[0][1]+k[1][0])*slope*intercept + k[1][1]*intercept**2
    constant_sign(obstruction)
    # 只保持循环置换时，六点分成两个三点轨道，各有同一实二次方程。
    rotation = [points.index((p[1], p[2], p[0])) for p in points]
    assert {4, rotation[4], rotation[rotation[4]]} == {4, 12, 13}
    assert {15, rotation[15], rotation[rotation[15]]} == {5, 10, 15}
    swap = [points.index((p[2], p[1], p[0])) for p in points]
    assert swap[4] == 15
    assert k[1][0] == k[0][1]
    constant_sign(h[1][1])
    # 不同根的和与积由韦达公式固定，因此 (4,15) 的内积仍给出矛盾。
    mixed = k[0][0]*h[1][1] - k[0][1]*(h[0][1]+h[1][0]) + k[1][1]*h[0][0]
    constant_sign(mixed)
    print("C3-invariant rank five excluded for every cutoff below 392.")
    print("Cyclic symmetry on the six outer points alone excludes arbitrary complex inner coordinates.")


if __name__ == "__main__":
    verify()
