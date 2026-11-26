"""前527枚实际零点的平方读数，加上归一平方，可线性恢复全部轨道系数乘积。"""
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    _, _, rows = run_path(str(root / "ch39-proton-detection-prefix.py"))["samples"](root, 527)
    pairs = list(combinations_with_replacement(range(32), 2))
    assert len(pairs) == len(rows)+1
    system = arb_mat([[arb(1)] * len(pairs)]
                     + [[row[i]*row[j] for i, j in pairs] for row in rows])
    determinant = system.det()
    assert not determinant.contains(0)
    print("归一平方与前527枚的528阶乘积矩阵行列式",determinant.str(12))


if __name__ == "__main__":
    verify()
