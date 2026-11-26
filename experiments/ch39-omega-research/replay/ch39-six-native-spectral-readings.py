"""检验原壳谱在六个高度上的唯一性，以及任选五项时的一维歧义。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    geometry = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 320
    heights = [zero.imag for zero in acb.zeta_zeros(2_306_888_100_085, 8)]
    offsets = [arb(0)] + [height - heights[0] for height in heights[1:]]
    assert all(a < b for a, b in zip(offsets, offsets[1:]))
    assert offsets[-1] < arb("1.61")
    angles = [(1 - arb(d) / 1734).acos() for d in (244, 724, 804)]
    frequencies = [tip - angles[0] for tip in angles[1:]]
    w1, w2 = frequencies
    initial = [w * heights[0] for w in frequencies]
    factors = [-36 * geometry["unit"](pair).imag for pair in ((510, 1), (1173, -2))]
    assert all(not factor.contains(0) for factor in factors)

    # 以有理网格覆盖整个高度区间，验证两个已知相位分量的商严格单调。
    second_lower = []
    wronskian_lower = []
    for k in range(161):
        offset = arb(k).union(arb(k + 1)) / 100
        a, b = [phase + w * offset for phase, w in zip(initial, frequencies)]
        second_lower.append(b.sin().lower())
        wronskian_lower.append((w1 * a.cos() * b.sin() - w2 * b.cos() * a.sin()).lower())
    assert min(second_lower) > arb(".138")
    assert min(wronskian_lower) > arb(".059")

    # 固定第二频率，截取四个振幅及第一频率的导数列。
    first_cos = factors[0] * initial[0].sin()
    first_sin = factors[0] * initial[0].cos()
    rows = []
    for x in offsets:
        a, b = w1 * x, w2 * x
        rows.append([a.cos(), a.sin(), b.cos(), b.sin(),
                     x * (-first_cos * a.sin() + first_sin * a.cos())])
    minors = [(abs(arb_mat([rows[i] for i in nodes]).det()), nodes)
              for nodes in combinations(range(8), 5)]
    assert len(minors) == 56 and all(det > 0 for det, _ in minors)
    print("全窗第二正弦、Wronskian 下界", min(second_lower).str(14), min(wronskian_lower).str(14))
    print("56 个五阶子式的绝对值下界", min(det.lower() for det, _ in minors).str(14))


if __name__ == "__main__":
    verify()
