"""从实际正零点的带符号质子读数恢复32个面轨道概率，并精确纠正一个错误读数。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    _, _, rows = run_path(str(root / "ch39-proton-detection-prefix.py"))["samples"](root, 46)
    ones = [arb(1)]*32
    # 取前31个带符号读数作基坐标，消去精确已知的归一行。
    inverse = arb_mat([ones] + rows[:31]).inv()
    reduced = (arb_mat(rows[31:33]) * inverse).tolist()
    assert len(reduced) == 2
    a, b = (row[1:] for row in reduced)
    assert min(abs(z).lower() for z in a+b) > 0
    minors = [a[i]*b[j] - a[j]*b[i] for i, j in combinations(range(31), 2)]
    minor_floor = min(abs(v).lower() for v in minors)
    assert minor_floor > 0
    # 前32枚去掉一项的各个31项读取方案；求对每项独立绝对误差的精确行和。
    norms = []
    for missing in range(32):
        keep = [ones] + [row for j, row in enumerate(rows[:32]) if j != missing]
        recover = arb_mat(keep).inv()
        norms.append(max(sum(abs(recover[i, j]) for j in range(1, 32)) for i in range(32)))
    best = min(range(32), key=lambda j: float(norms[j]))
    assert all(norms[j] > norms[best] for j in range(32) if j != best)
    assert best == 0 and norms[best] < arb("24085.843") and norms[-1] < arb("40460.906")
    # 前16枚可被严格正的轨道概率同时消零。其余30行加归一最多31条约束。
    support = (0, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 26)
    outside = sorted(set(range(32)) - set(support))
    normalization = [ones] + rows[:16]
    small = arb("1e-8")
    base = arb_mat([[row[j] for j in support] for row in normalization])
    rhs = arb_mat([[int(i == 0) - small*sum(row[j] for j in outside)]
                   for i, row in enumerate(normalization)])
    assert min(base.solve(rhs).entries()) > small
    print("前46枚绝对值碰撞见证中点的最小轨道概率",small.str(12))
    print("前33枚的双行子式绝对值下界", minor_floor.str(17))
    print("前32枚中略去的最佳编号、误差放大、原31枚前缀误差放大",
          best+1, norms[best].str(17), norms[-1].str(17))


if __name__ == "__main__":
    verify()
