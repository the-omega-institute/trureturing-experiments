"""面系数改为归一化实或复振幅后，32项实读出矩阵与自然计数Gram矩阵的锐界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def three_zero_cancellation(root):
    rows = run_path(str(root / "ch39-unrestricted-proton-detector.py"))["verify"](root)
    columns = [[arb(1)] + [row[j] for row in rows] for j in range(32)]
    # 非负部用轨道5、27、31，负部用轨道28。最小化负部质量。
    basis = [(5, 1), (27, 1), (31, 1), (28, -1)]
    matrix = arb_mat([[sign*columns[j][i] for j, sign in basis] for i in range(4)])
    primal = matrix.solve(arb_mat([[1], [0], [0], [0]])).entries()
    dual = matrix.transpose().solve(arb_mat([[0], [0], [0], [1]])).entries()
    assert min(primal) > 0 and (primal[-1] - dual[0]).contains(0)
    for j, column in enumerate(columns):
        price = sum(a*b for a, b in zip(column, dual))
        if (j, 1) not in basis:
            assert price < 0
        if (j, -1) not in basis:
            assert -price < 1
    assert arb(".03201774") < primal[-1] < arb(".03201775")
    print("使三个实际编号同时为零的最小负面系数总量", primal[-1].str(17))


def verify():
    root = Path(__file__).resolve().parent
    _, _, rows = run_path(str(root / "ch39-proton-detection-prefix.py"))["samples"](root, 32)
    matrix = arb_mat(rows)
    ones = arb_mat([[1]] * 32)
    dual = matrix.transpose().solve(ones).entries()
    gap = min(abs(v).lower() for v in dual)
    assert gap > 0
    value = 1/sum(abs(v) for v in dual)
    signs = arb_mat([[1 if v > 0 else -1] for v in dual])
    extremal = matrix.solve(signs) * value
    assert (sum(extremal.entries())-1).contains(0) and arb(".00230079323") < value < arb(".00230079324")
    print("32个编号的反向权绝对值下界、可达最大模下界", gap.str(15), value.str(17))
    three_zero_cancellation(root)
    gram = run_path(str(root / "ch39-proton-all-zero-rms.py"))["verify"](root)
    sol = gram.solve(ones)
    variance = 1/sum(sol.entries())
    optimum = sol*variance
    negative_mass = sum(-v for v in optimum.entries() if v < 0)
    assert all(v > 0 or v < 0 for v in optimum.entries()) and negative_mass > 0
    assert arb(".01055307094") < variance.sqrt() < arb(".01055307095")
    print("归一复面系数的全谱最小均方根、最优解负系数绝对值和",
          variance.sqrt().str(17), negative_mass.str(17))


if __name__ == "__main__":
    verify()
