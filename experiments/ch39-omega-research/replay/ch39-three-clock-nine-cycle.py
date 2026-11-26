"""用九环的二次代数数解实现截止 392 的最小辅助维数。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat, fmpq_poly

pairs = run_path(str(Path(__file__).with_name("ch39-three-clock-outer-pairs.py")))
previous = pairs["previous"]
submatrix = pairs["submatrix"]


def verify():
    _, _, _, g0, g1, g2 = previous["MEM"]["coefficients"]()
    delta = g2 - g1 * g0.solve(g1)
    outer, centre = pairs["outer"], pairs["centre"]
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17, 17, 17)})
    n = len(points)
    def distance(i, j):
        return local["squared"](points[i], points[j])

    t = fmpq(1, 10**10)
    gram = pairs["projected_outer_gram"](delta, t)
    metric = fmpq_mat(7, 7)
    for i in range(6):
        for j in range(6):
            metric[i, j] = gram[i, j]
    metric[6, 6] = t

    # 第 i 行向量写成 a_i + y_i b_i。外点与中点固定。
    a, b = fmpq_mat(n, 7), fmpq_mat(n, 7)
    for k, i in enumerate(outer):
        a[i, k] = a[i, 6] = 1
    a[centre, 6] = delta[centre, outer[0]] / t
    inner = sorted(set(range(n)) - set(outer) - {centre})
    for i in inner:
        fixed = [k for k, j in enumerate(outer) if distance(i, j) > 392]
        assert len(fixed) == 4
        rhs = fmpq_mat([[delta[i, outer[k]], -t] for k in fixed])
        cols = submatrix(gram, fixed, fixed).solve(rhs)
        for l, k in enumerate(fixed):
            a[i, k], b[i, k] = cols[l, 0], cols[l, 1]
        b[i, 6] = 1

    aa, ba, bb = a * metric * a.transpose(), b * metric * a.transpose(), b * metric * b.transpose()
    x = fmpq_poly([0, 1])
    orbit_x = {4, 5, 10, 12, 13, 15}
    orbit_z = {6, 7, 16}
    # 第 4、15 点的约束给出 P(x)=0。第 4、16 点的约束给出 z=u(x)/v(x)。
    p = (aa[4, 15] - delta[4, 15]) + (ba[4, 15] + ba[15, 4]) * x + bb[4, 15] * x**2
    u = delta[4, 16] - aa[4, 16] - ba[4, 16] * x
    v = ba[16, 4] + bb[4, 16] * x
    left, right = fmpq(-6302, 1000), fmpq(-6301, 1000)
    assert p.degree() == 2 and p(left) * p(right) < 0
    assert v(left) * v(right) > 0
    lower_den = min(abs(v(left)), abs(v(right)))**2
    y_num = [x * v if i in orbit_x else u if i in orbit_z else fmpq_poly()
             for i in range(n)]
    den = v**2
    row_bounds = [fmpq(0) for _ in range(n)]
    residues = fmpq_mat(n, n), fmpq_mat(n, n)
    for i in range(n):
        for j in range(i, n):
            num = aa[i, j] * den + (ba[i, j] * y_num[i] + ba[j, i] * y_num[j]) * v
            num += bb[i, j] * y_num[i] * y_num[j]
            if distance(i, j) > 392:
                assert (num - delta[i, j] * den) % p == 0
            line = num % p
            for degree, m in enumerate(residues):
                m[i, j] = m[j, i] = line[degree]
            bound = max(abs(line(left)), abs(line(right))) / lower_den
            row_bounds[i] += bound
            if i != j:
                row_bounds[j] += bound
    for perm in permutations(range(3)):
        idx = [points.index(tuple(p[k] for k in perm)) for p in points]
        assert all(m == submatrix(m, idx, idx) for m in residues)
    norm = max(row_bounds)
    assert norm < fmpq(1, 10**7)
    inverse_norm = max(sum(abs(z) for z in row) for row in g0.inv().tolist())
    assert inverse_norm * norm < (1 - fmpq(98, 1000))**2
    print("cutoff 392, auxiliary dimension 24, certified row norm below", float(norm))


if __name__ == "__main__":
    verify()
