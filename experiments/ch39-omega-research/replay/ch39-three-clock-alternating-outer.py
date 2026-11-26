"""六外点交替方向首次可落入三时钟 Gram 核的精确平方弦距截止。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

plane = run_path(str(Path(__file__).with_name("ch39-three-clock-plane-completion.py")))


def verify():
    points, distance, _, delta, _, _ = plane["shell_data"]()
    group = list(permutations(range(3)))
    signs = [(-1)**sum(g[a] > g[b] for a, b in combinations(range(3), 2))
             for g in group]

    def orbit(rep):
        return [points.index(tuple(points[rep][k] for k in g)) for g in group]

    outside, inside = orbit(0), orbit(4)
    assert len(set(outside)) == len(set(inside)) == 6
    block = plane["pairs"]["submatrix"](delta, inside, outside)
    sign_column = fmpq_mat([[s] for s in signs])
    eigenvalue = (block*sign_column)[0, 0]
    assert block*sign_column == eigenvalue*sign_column
    assert fmpq(6497739194, 10**10) < eigenvalue < fmpq(6497739195, 10**10)
    for k, i in enumerate(inside):
        assert distance[i][outside[k]] == 216
        assert all(distance[i][j] >= 296 for j in outside if j != outside[k])

    # 截止 216 的构造在外点上使用 t=1 的投影 Gram 加全一矩阵。
    outer = plane["outer"]
    centered = plane["pairs"]["projected_outer_gram"](delta, 1)
    positions = [outer.index(i) for i in outside]
    gram = fmpq_mat([[centered[a, b] + 1 for b in positions] for a in positions])
    assert gram*sign_column == fmpq_mat(6, 1)
    return points, distance, outside, inside, sign_column, eigenvalue


if __name__ == "__main__":
    *_, eigenvalue = verify()
    print("First cutoff for a null alternating outer direction: 216.",
          "Scaled alternating coupling:", float(eigenvalue))
