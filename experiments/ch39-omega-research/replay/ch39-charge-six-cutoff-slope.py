"""荷六第一枢轴要求在很窄的两个原壳能级之间拉开能标。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    fields = run_path(str(root / "ch39-three-index-free-atoms.py"))
    _, heat, (rows,) = fields["broken_fields"](root, [2_306_888_100_085])
    ctx.prec = 256
    _, levels, direct, _ = next(item for item in native["shell_weights"](root) if item[0] == 6)
    weights = [v.real for v in direct[0]]
    tv = sum(abs(c) for c in weights)
    alpha = -sum(a for _, a, _ in rows) / (2 * arb.pi())
    mass_error = tv * sum(abs(a) * (1 - heat(2, z)) for _, a, z in rows) / (2 * arb.pi())
    neg = sum(max(arb(0), -c.lower()) for c in weights[:-1])
    assert weights[-1] > 0 and mass_error < 1
    x18 = (1 - mass_error) / (alpha * neg)
    # 所有负权的单节点弱读数均由最大许可位置的下界控制。
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    floor = native["weak_series"]()(arb(1))
    early = sum(max(arb(0), c.upper()) for c in weights[:15])
    late = sum(max(arb(0), c.upper()) for c in weights[15:])
    noise = tv * (arb("2.398e-21") + arb("7.468e-36"))
    room = -ref - early * (floor - ref) - noise
    assert room > 0 and neg > 0
    x16 = late * arb("7.468e-36") / room
    slope = (x18/x16).log() / (levels[17]/levels[15]).log()
    print("第16层平方截止上界、第18层下界、本征值比、必要的双对数割线斜率",
          x16.upper().str(15), x18.lower().str(15),
          (levels[17]/levels[15]).upper().str(15), slope.lower().str(15))
    assert x18 > arb(".14330") and x16 < arb("4.26148e-34")
    assert slope > arb("4626.83") and levels[17]/levels[15] < arb("1.01632")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
