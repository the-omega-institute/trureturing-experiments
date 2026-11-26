"""前32枚实际零点里只随机测一枚时，最优的保证均方根及唯一抽样律。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    _, _, rows = run_path(str(root / "ch39-proton-detection-prefix.py"))["samples"](root, 32)
    # 正负编号表示要将该行乘以 +1 或 -1。轨道列使用从零开始的编号。
    signed_indices = (4, 14, 16, 20, 23, 27, -5, -7, -13, -15, -17, -22, -24, -25, -31)
    support = (2, 6, 7, 11, 12, 13, 14, 17, 18, 19, 20, 23, 25, 26, 29)
    active = {abs(n) for n in signed_indices}
    assert len(active) == len(support) == len(signed_indices)
    signs = [1 if n > 0 else -1 for n in signed_indices]
    signed = [[s * v for v in rows[abs(n)-1]] for s, n in zip(signs, signed_indices)]
    # 令每个活动观测的有符号值等于 z，并将面概率归一。
    table = arb_mat([[row[j] for j in support] + [-1] for row in signed]
                    + [[1] * len(support) + [0]])
    primal = table.solve(arb_mat([[0]] * len(support) + [[1]])).entries()
    masses, level = primal[:-1], primal[-1]
    dual = table.transpose().solve(arb_mat([[0]] * len(support) + [[-1]])).entries()
    probabilities = dual[:-1]
    assert min(masses + probabilities) > 0
    assert (sum(probabilities) - 1).contains(0) and (level + dual[-1]).contains(0)
    assert arb(".0108810990718") < level < arb(".0108810990719")
    for n, row in enumerate(rows, 1):
        if n not in active:
            assert abs(sum(p * row[j] for j, p in zip(support, masses))) < level
    for j in range(32):
        if j not in support:
            assert sum(prob * row[j] for row, prob in zip(signed, probabilities)) > level
    print("最优最坏均方根、平均绝对值及其平方", level.str(18), (level**2).str(18))
    for n, prob in sorted(zip(signed_indices, probabilities), key=lambda item: abs(item[0])):
        print("带符号编号、最优抽样概率", n, prob.str(14))


if __name__ == "__main__":
    verify()
