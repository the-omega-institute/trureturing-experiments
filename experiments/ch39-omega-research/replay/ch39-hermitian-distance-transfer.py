"""十九站原生弦距传输取自伴系数时的最佳条件重泄漏。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def real_gram(transport, outputs):
    columns = transport["columns"](outputs)
    gram = transport["adjoint"](columns) * columns
    return acb_mat([[gram[i, j].real for j in range(gram.ncols())]
                    for i in range(gram.nrows())])


def verify():
    root = Path(__file__).resolve().parent
    transport = run_path(str(root / "ch39-joint-zero-distance-class-transfer.py"))
    positive_definite = run_path(str(root / "ch39-joint-zero-diagonal-light-filters.py"))["certify_positive"]
    light, _, _, height, center, identity, _, outer, spokes, source, _ = transport["shell_classes"](root)
    classes = outer + list(spokes.values())
    assert len(classes) == 24
    outputs = [matrix * source for matrix in classes]
    total = real_gram(transport, outputs)
    positive_definite(total)
    for charge, lower, upper in ((2, ".0000589009", ".0000589011"),
                                 (-4, ".0001713", ".0001714")):
        target = light["light_basis"](charge, height, center)
        heavy = identity - target * transport["adjoint"](target)
        loss = real_gram(transport, [heavy * matrix for matrix in outputs])
        lo, hi = arb(lower), arb(upper)
        positive_definite(loss - total * lo)
        assert (loss - total * hi).det().real < 0
        print("目标荷、实弦距系数的最小条件重泄漏区间", charge, lower, upper)


if __name__ == "__main__":
    verify()
