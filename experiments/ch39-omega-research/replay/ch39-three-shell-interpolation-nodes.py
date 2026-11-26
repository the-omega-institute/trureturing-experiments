"""核验三枚相邻零点的三荷谱，共同有界插值所需的 171 个节点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    ctx.prec = 320
    center = 2_306_888_100_085
    nodes = []
    for number in (center - 1, center, center + 1):
        height = acb.zeta_zero(number).imag
        for charge in (-4, 2, 6):
            levels = ordered(shell, charge, height)
            assert len(levels) == 19
            nodes.extend((level, number, charge, i) for i, level in enumerate(levels, start=1))
    nodes.sort(key=lambda row: row[0].mid())
    gaps = [(b[0] - a[0], a[1:], b[1:]) for a, b in zip(nodes, nodes[1:])]
    assert len(nodes) == 171 and all(gap > arb("8.28787e-7") for gap, _, _ in gaps)
    gap, left, right = min(gaps, key=lambda row: row[0].mid())
    assert all(gap < other for other, a, b in gaps if (a, b) != (left, right))
    print("节点数、最小间距与两端的（编号、荷、层）", len(nodes), gap.str(16), left, right)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
