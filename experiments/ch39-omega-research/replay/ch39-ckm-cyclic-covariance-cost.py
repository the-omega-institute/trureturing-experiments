"""参考 CKM 矩阵要求的三阶轮换协变性偏差，容许质量线独立改名和调相。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    center = ckm["standard_matrix"](*ckm["centers"])
    window = ckm["target_box"]()
    for matrix, name in ((center, "参考中心"), (window, "原四参数容差盒")):
        # 允许源、目标各自重排；迹中只出现三对相同特征标的相对匹配。
        kept = [sum(abs(matrix[i][j])**2 for i,j in enumerate(perm))
                for perm in permutations(range(3))]
        assert all(kept[0] > v for v in kept[1:])
        off = sum(abs(matrix[i][j])**2 for i in range(3) for j in range(3) if i!=j)
        error = (3*off).sqrt()
        print(name,"最佳相同特征标匹配、单映射 Frobenius 偏差、分摊时至少一支",
              kept[0].str(15),error.str(15),(error/2).str(15))
        assert error > arb(".560")


if __name__ == "__main__":
    verify()
