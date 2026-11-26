"""五壳三荷 285 处互异谱输入及一个共同有界实有理响应的存在条件。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root):
    center = 2_306_888_100_085
    floor = arb("91.1876") / arb("1.22089e19")
    _, flavor = run_path(str(root / "ch39-five-shell-three-space-quark-witness.py"))["verify"](root)
    inner, _ = run_path(str(root / "ch39-three-shell-two-space-readout.py"))["verify"](root, return_data=True)
    outer = run_path(str(root / "ch39-outer-monotone-distance-gravity.py"))["verify"](root, training=inner)
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    nodes = []
    for n in range(center - 2, center + 3):
        ctx.prec = 320
        height = acb.zeta_zero(n).imag
        for q in (-4, 2, 6):
            levels = ordered(shell, q, height)
            if q != 6:
                values = flavor[n, q]
            elif n in outer:
                values = outer[n]
            else:
                values = []
                for level in levels:
                    matches = [y for energy, y in inner if level.overlaps(energy)]
                    assert len(matches) == 1
                    values.extend(matches)
            assert len(values) == len(levels) == 19
            assert all(floor < y < 1 for y in values)
            nodes.extend((energy, n, q, i) for i, energy in enumerate(levels, start=1))
    nodes.sort(key=lambda row: row[0].mid())
    assert len(nodes) == 285
    gaps = [(b[0] - a[0], a[1:], b[1:]) for a, b in zip(nodes, nodes[1:])]
    assert all(gap > 0 for gap, _, _ in gaps)
    gap, left, right = min(gaps, key=lambda row: row[0].mid())
    assert all(gap < v for v, a, b in gaps if (a, b) != (left, right))
    print("五壳谱点数、最小间距及对应的编号荷与层", len(nodes), gap.str(13), left, right)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
