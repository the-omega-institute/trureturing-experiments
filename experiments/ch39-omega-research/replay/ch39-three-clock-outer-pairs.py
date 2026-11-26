"""十九点壳三时钟的两个精确辅助维数。只用有理数与有理多项式。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat, fmpq_poly

previous = run_path(str(Path(__file__).with_name("ch39-local-memory-rank.py")))
submatrix = previous["submatrix"]
outer = (0, 2, 1, 3, 17, 18)
centre = 11


def det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def adj2(m):
    return [[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]]


def multiply2(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def positive_on_ray(p):
    assert p[0] > 0 and all(c >= 0 for c in p.coeffs())


def rank_obstruction(delta):
    # 投掉中点后，三个跨块的每个条目都减去同一 t。无中点时取 t=0。
    a, b, c = outer[:2], outer[2:4], outer[4:]
    def block(rows, columns):
        return [[fmpq_poly([delta[i, j], -1]) for j in columns] for i in rows]

    ba, bc, ac = block(b, a), block(b, c), block(a, c)
    den = det2(ba) * det2(ac)
    num = multiply2(multiply2(adj2(ba), bc), adj2(ac))
    # -2 den² Herm(BA^{-1} BC AC^{-1}) 的实对称矩阵。
    h = [[-den * (num[i][j] + num[j][i]) for j in range(2)] for i in range(2)]
    positive_on_ray(h[0][0])
    positive_on_ray(det2(h))


def projected_outer_gram(delta, offset):
    signs = (1, -1, -1, 1, 1, -1)
    gram = submatrix(delta, outer, outer)
    for i in range(6):
        for j in range(6):
            if i // 2 != j // 2:
                gram[i, j] -= offset
    for i in range(6):
        other = [j for j in range(6) if i // 2 != j // 2]
        gram[i, i] = -sum(gram[i, j] for j in other if signs[i] == signs[j])
        gram[i, i ^ 1] = -sum(gram[i, j] for j in other if signs[i] != signs[j])
    assert gram == gram.transpose() and gram.rank() == 4
    assert all(submatrix(gram, range(k), range(k)).det() > 0 for k in range(1, 5))
    return gram


def completion(delta, distances, cutoff):
    n = delta.nrows()
    offset = fmpq(1, 10**10) if cutoff < 544 else fmpq(0)
    gram = projected_outer_gram(delta, offset)
    w = fmpq_mat(n, n)
    for k, i in enumerate(outer):
        for l, j in enumerate(outer):
            w[i, j] = gram[k, l] + offset
    if offset:
        for i in outer:
            w[i, centre] = w[centre, i] = delta[i, centre]
        w[centre, centre] = delta[centre, outer[0]]**2 / offset

    # 依次在已经构造的向量中，取满足全部已知内积的最小范数向量。
    done = [*outer, centre]
    for i in sorted(set(range(n)) - set(done)):
        fixed = [j for j in done if distances[i, j] > cutoff]
        if fixed:
            block = submatrix(w, fixed, fixed)
            rhs = submatrix(delta, fixed, [i])
            coeff = block.solve(rhs)
            row = coeff.transpose() * submatrix(w, fixed, done)
            for k, j in enumerate(done):
                w[i, j] = w[j, i] = row[0, k]
            w[i, i] = (rhs.transpose() * coeff)[0, 0]
        done.append(i)
    assert all(w[i, j] == delta[i, j] for i, j in combinations(range(n), 2)
               if distances[i, j] > cutoff)
    assert w.rank() == (5 if offset else 4)
    return w


def verify():
    _, _, _, g0, g1, g2 = previous["MEM"]["coefficients"]()
    delta = g2 - g1 * g0.solve(g1)
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17, 17, 17)})
    distances = {(i, j): local["squared"](p, q)
                 for i, p in enumerate(points) for j, q in enumerate(points)}
    assert min(distances[i, j] for k, i in enumerate(outer)
               for l, j in enumerate(outer) if k // 2 != l // 2) == 1152
    assert {distances[centre, i] for i in outer} == {544}
    assert len({delta[centre, i] for i in outer}) == 1 and delta[centre, outer[0]] != 0
    rank_obstruction(delta)
    inverse_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    assert inverse_norm < 315000
    for cutoff in (504, 512, 544):
        w = completion(delta, distances, cutoff)
        norm = max(sum(abs(x) for x in row) for row in w.tolist())
        assert inverse_norm * norm < (1 - fmpq(98, 1000))**2
        print("cutoff", cutoff, "rank", w.rank(), "auxiliary states", 19 + w.rank(),
              "absolute row sum", float(norm))


if __name__ == "__main__":
    verify()
