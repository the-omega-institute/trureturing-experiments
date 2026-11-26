"""第 2 与第 30 枚真实零点共同引力读数所需的最少有效截止翻号次数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def verify(root):
    # 零次翻号的下界用两个真实背景的全正半轴严格导数序。
    run_path(str(root / "ch39-two-thirty-universal-cutoff.py"))["verify"]()
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, spectra = fields(root, (2, 30))

    def moments(x):
        (k, w), (other_k, other_w) = [single["moments"](rows, lib["expint"], x) for rows in spectra]
        return [k, w, other_k-k, other_w-w]

    ends = list(map(arb, ("6.5e-9", "1.2e-8", "1.35e-8", "1.55e-8", ".01")))
    readings = [moments(x) for x in ends]
    target = [1, arb("1e-122"), 0, 0]
    # H 在最初一段恒为 1。矩阵的列依次为其余四段的热核导数积分。
    matrix = arb_mat([[b-a for a, b in zip(left, right, strict=True)]
                      for left, right in zip(readings[:-1], readings[1:], strict=True)]).transpose()
    rhs = arb_mat([[t-a] for t, a in zip(target, readings[0], strict=True)])
    det = matrix.det()
    assert arb("2.6e-41") < det < arb("2.8e-41")

    def midpoints(m):
        return arb_mat([[m[i, j].mid() for j in range(m.ncols())] for i in range(m.nrows())])

    pre = midpoints(midpoints(matrix).inv())
    assert not pre.det().contains(0)
    initial = midpoints(pre*rhs)
    heights = (initial + (pre*matrix).solve(pre*(rhs-matrix*initial))).entries()
    boxes = (("4.02e7", "4.04e7"), ("-4.12e7", "-4.10e7"),
             ("-2.15e6", "-2.05e6"), ("-3.272e-5", "-3.270e-5"))
    assert all(arb(lo) < h < arb(hi) for (lo, hi), h in zip(boxes, heights, strict=True))

    # 相邻段高度之差就是截止测度的原子权，望远镜求和等于 1。
    profile = [arb(1), *heights, arb(0)]
    weights = [a-b for a, b in zip(profile[:-1], profile[1:], strict=True)]
    assert sum(weights).contains(1)
    assert all(sum(p*row[i] for p, row in zip(weights, readings, strict=True)).contains(t)
               for i, t in enumerate(target))
    print("四段系统的行列式与剖面高度", det.str(15), *[h.str(15) for h in heights])
    print("五个截止原子的权重", *[p.str(15) for p in weights])
    print("最少翻号次数为 1；翻号处的 x 与最大 Q/M_P", ends[1].str(12), ends[-1].sqrt().str(12))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
