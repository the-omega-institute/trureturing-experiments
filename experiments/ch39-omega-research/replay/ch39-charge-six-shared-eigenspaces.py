"""保持七维不变空间的荷六输运，其立方根跨空间转移的共同门槛与两个锐例。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    transfer = run_path(str(root / "ch39-quark-cyclic-subspace-transfer.py"))
    transfer["verify"](root)
    ctx.prec = 320
    threshold = arb(".50134")
    omega = acb(-1, arb(3).sqrt()) / 2
    # 一条内线、一条外线共享荷六本征值。相差 120 度的两条根线作 45 度混合。
    split = acb_mat([[1, 1], [1, -1]]) * (1 / arb(2).sqrt())
    root2 = split * acb_mat([[1, 0], [0, omega]]) * adj(split)
    identity2 = acb_mat([[1, 0], [0, 1]])
    assert (root2**3 - identity2).contains(acb_mat(2, 2))
    assert ((root2**2).conjugate().transpose() - root2).contains(acb_mat(2, 2))
    probability = abs(root2[1, 0])**2
    assert probability.contains(arb(3) / 4) and probability > threshold
    # 混合本征空间至少三维时，选跨内外的一个三循环。
    cycle = acb_mat([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    identity3 = acb_mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert cycle**3 == identity3 and cycle[1, 0] == 1
    # 简单谱情形下违反对称所需的 Frobenius 交换子下界。
    commutator = (threshold / 2).sqrt()
    assert commutator > arb(".50066")
    print("两荷共同概率门槛、二维与高维共享本征空间最大值、简单谱交换子比下界",
          threshold.str(8), probability.str(10), 1, commutator.str(12), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
