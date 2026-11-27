"""固定谱底时，七壳弱差与真空项对任意秩共享正源的超额迹上界。"""

from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def trace_bound(shells, floor, ids, caps):
    planck = arb("1.22089e19")
    radii = [(arb(q) / planck - floor).sqrt().upper() for q in caps]
    selected = acb_mat([shell[0].tolist()[i] for shell, i in zip(shells, ids, strict=True)])
    inverse = selected.inv()
    column_norms = [
        sum(abs(inverse[k, j]) ** 2 for k in range(7)).sqrt() for j in range(7)
    ]
    trivial_trace = sum(v * r for v, r in zip(column_norms, radii, strict=True)) ** 2

    thresholds = []
    for _, pos, positive_rows, weights, c, error, allowance, holder in shells:
        positive = arb(0)
        for i, row in zip(pos, (positive_rows * inverse).tolist(), strict=True):
            row_norm = sum(abs(z) * r for z, r in zip(row, radii, strict=True))
            # 一般正源的最大本征值可以任意大，此处只用七条低线的约束。
            amplitude = (floor + min(trivial_trace.upper(), (row_norm**2).upper())).upper()
            positive += weights[i] * (c * amplitude**4 + error * amplitude**2)
        cost = arb(".99") * c * holder
        assert cost > 0 and positive + allowance > 0
        # 余下十二条谱行的截止和等于 12f 加上该子空间的超额迹。
        complementary_cap = ((positive + allowance) / cost).sqrt().sqrt()
        thresholds.append((trivial_trace + complementary_cap - 12 * floor).upper())
    return min(thresholds)


def experiment(root, levels):
    prepare = run_path(str(root / "ch39-seven-unrestricted-vacuum-weak.py"))["prepare"]
    candidates, shells, refined = prepare(root, levels)
    floor, heavy, top = levels
    bounds = []
    for ids in product(*candidates):
        caps = refined if ids == (0,) * 7 else (3000,) * 7
        bounds.append((trace_bound(shells, floor, ids, caps), ids))
    assert len(bounds) == 972
    limit, worst = max(bounds, key=lambda item: item[0])
    original = heavy + top - 2 * floor
    deficit = original - limit
    relative = deficit / original
    assert worst == (0,) * 7
    assert limit < arb("1.13225") and original > arb("1.185088")
    assert deficit > arb(".052838") and relative > arb(".04458")
    print("七壳选线数、最大迹上界及其一基选线", len(bounds), limit.str(16), tuple(i + 1 for i in worst))
    print("原超额迹、必要赤字下界、相对赤字下界", original.str(16),
          deficit.lower().str(16), relative.lower().str(16), flush=True)
    return limit


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    levels = source["verify"](root, return_operator="spectral")[2]
    experiment(root, levels)
