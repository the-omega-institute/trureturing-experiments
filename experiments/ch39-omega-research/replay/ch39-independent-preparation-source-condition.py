"""独立改变上下型制备时，夸克质量窗对三源换基条件数的全域下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    engine = run_path(str(Path(__file__).with_name("ch39-shared-724-quark-gaps.py")))
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    # 上型有效第一角为 -2x；x 属于 [0,π] 已绕整圆一周。
    count = engine["verify"](-4, 1, 0, arb.pi(), 1, "8.1", source_basis=identity)
    target = arb(".8573")/arb(".00316")
    lower = target/(arb("8.1")*9/7)
    assert lower > arb("26.05")
    print("上型全域终格数、正交化三源坐标中换基条件数的下界", count, lower.str(12))


if __name__ == "__main__":
    verify()
