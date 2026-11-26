"""既定十七位置、十七权与门户对四个留出实际零点的引力读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    witness = run_path(str(root / "ch39-three-index-positive-running.py"))
    _, weights, _, _ = witness["corrected"](root)
    sites = list(map(arb, witness["SITES"]))
    center = 2_306_888_100_085
    offsets = (0, -2, 2, -3, 3)
    broken_fields = run_path(str(root / "ch39-three-index-free-atoms.py"))["broken_fields"]
    moments, heat, fields = broken_fields(root, [center + j for j in offsets])
    # 在相乘前逐列扣除已匹配的中心矩，保留公共质量谱相消带来的精度。
    baseline = [moments(fields[0], heat, x) for x in sites]
    for offset, rows in zip(offsets[1:], fields[1:], strict=True):
        measured = [moments(rows, heat, x) for x in sites]
        newton, vacuum = [
            sum(c * (v[i] - v0[i]) for c, v, v0 in zip(weights, measured, baseline, strict=True))
            for i in (0, 1)
        ]
        print("零点编号、偏移、K−1、W−目标", center + offset, offset,
              newton.str(17), vacuum.str(17), flush=True)
        assert newton < arb("-0.00043710") and vacuum > arb("7.88e-10")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
