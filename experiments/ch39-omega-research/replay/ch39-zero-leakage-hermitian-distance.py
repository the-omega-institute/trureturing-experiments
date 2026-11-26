"""三代零重泄漏传输距共同相位自伴矩阵的最小相对 Frobenius 距离。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


# 各目标重投影下用作基的十六列，按自环、二十种非辐边、三种辐边排序。
CASES = (
    (2, (17, 5, 1, 6, 13, 8, 18, 16, 4, 14, 10, 7, 15, 23, 19, 12),
     ".00000043130", ".00000043135"),
    (-4, (17, 5, 1, 6, 9, 8, 18, 4, 16, 14, 10, 15, 7, 23, 13, 11),
     ".0000013880", ".0000013882"),
)


def verify():
    root = Path(__file__).resolve().parent
    transport = run_path(str(root / "ch39-joint-zero-distance-class-transfer.py"))
    positive = run_path(str(root / "ch39-joint-zero-diagonal-light-filters.py"))["certify_positive"]
    adjoint = transport["adjoint"]
    columns = transport["columns"]
    light, cyclic, symmetry, height, center, eye, _, main, spokes, source, _ = transport["shell_classes"](root)
    classes = main + [spokes[d] for d in sorted(spokes)]
    # 支撑互斥的 24 个自伴类矩阵作 Frobenius 归一。虚系数的范数即反自伴部分的范数。
    outputs = [(matrix * source) / sum(abs(z)**2 for z in matrix.entries()).sqrt()
               for matrix in classes]
    identity = acb_mat([[int(i == j) for j in range(24)] for i in range(24)])
    first_coordinate = acb_mat([[int(i == 0)] for i in range(24)])
    for charge, selected, lower, upper in CASES:
        target = light["light_basis"](charge, height, center)
        assert sorted(cyclic["character_order"](target, symmetry)) == [0, 1, 2]
        heavy = eye - target * adjoint(target)
        blocks = [heavy * a for a in outputs]
        rows = adjoint(columns([blocks[i] for i in selected])) * columns(blocks)
        metric = rows * adjoint(rows)
        positive(metric)
        # 行空间维数为十六。上界来自目标重空间的三阶轮换交织空间维数。
        row_projector = adjoint(rows) * metric.solve(rows)
        real_part = acb_mat([[row_projector[i, j].real for j in range(24)] for i in range(24)])
        lo, hi = arb(lower), arb(upper)
        shifted = real_part - identity * lo
        positive(shifted)
        # 平移反迭代一次，精确投影回零漏核，在给定总相位上直接检验距离上界。
        trial = (identity - row_projector) * shifted.solve(first_coordinate)
        norm_squared = sum(abs(z)**2 for z in trial.entries())
        fraction = sum(z.imag**2 for z in trial.entries()) / norm_squared
        assert fraction < hi
        out = sum((c * block for c, block in zip(trial.entries(), outputs, strict=True)),
                  acb_mat(19, 3))
        positive(adjoint(out) * out)
        print("目标荷、最小相对平方距离的下界与秩三见证上界",
              charge, lo.str(12), fraction.str(14))


if __name__ == "__main__":
    verify()
