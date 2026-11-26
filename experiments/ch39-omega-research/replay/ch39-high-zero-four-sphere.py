"""核验零软化 Weyl 点处第二份原壳的谱计数与共同能量余量。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    cone = run_path(str(root / "ch39-first-chern-edge.py"))
    sampling = run_path(str(root / "ch39-actual-zero-weyl-mesh.py"))
    levels = run_path(str(root / "ch39-weyl-light-yukawa-determinant.py"))["spectator_levels"]
    ctx.prec = 256
    height, energy, weights = sampling["flat_certificate"](native, cone)
    delta = arb(23)/25 * arb(2)**-20 / (height**2 + native["soft"]["point"]["study"]["offset"])
    assert len(native["soft"]["star"]["orbit"]) == 3
    common_lower, common_upper = [], []
    # 零软化在三条单位辐边上各加 delta，正半定改动的范数至多 6 delta。
    for charge, count in ((6, 3), (-4, 3), (2, 4), (1, 4)):
        eig = levels(native, cone["EDGES"], height, weights, charge)
        assert len(eig) == 19 and all(z+6*delta < energy for z in eig[:count])
        assert all(z > energy for z in eig[count:])
        lower_gap = energy - eig[count-1] - 6*delta
        upper_gap = eig[count] - energy
        positive_shift = arb(0).union(6*delta)
        common_lower.append(eig[count-1] + positive_shift)
        common_upper.append(eig[count] + positive_shift)
        assert lower_gap > arb(".05") and upper_gap > arb(".05")
        print("荷、低于交点能量的谱数、上下两侧能量余量下界",
              charge, count, lower_gap.lower().str(12), upper_gap.lower().str(12))
    # 分母至多 6 时，9 以下最大的分数是 53/6。
    first = arb(62)/7
    assert all(v < first for v in common_lower) and all(first < v for v in common_upper)
    assert any(v > arb(53)/6 for v in common_lower) and any(v < 9 for v in common_upper)
    print("共同可用的最小分母半径", "62/7", max(common_lower, key=float).str(10),
          min(common_upper, key=float).str(10))


if __name__ == "__main__":
    verify()
