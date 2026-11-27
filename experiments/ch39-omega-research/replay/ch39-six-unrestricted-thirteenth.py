"""冻结的六壳十八式见证在第十三壳原生第三枢轴上的读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def readout(root, levels, weighted):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    number = 2_306_888_100_098
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    adj = lib("common-three-shell-nine-endpoints")["adj"]
    gram = native["gram_cholesky"](root)
    _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
    weights = [
        z.real for z in native["readout"](basis, range(19), [14, 8, 7, 6], gram[3], 3)
    ]
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [number])

    floor = levels[0]
    rows = adj(basis) * weighted
    # 两列已吸收谱隙平方根，引力计算中按单位增量代入。
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    values, _ = gravity(rows, [], (floor, floor + 1, floor + 1), weights, heat, field, True)
    newton = values[0, 0] + 1
    vacuum = values[1, 0] + arb("1e-122")

    planck, mz = arb("1.22089e19"), arb("91.1876")
    p0, p1 = (mz / planck)**2, (arb("1.176397e14") / planck)**2
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    minimum = min((floor + sum(abs(z)**2 for z in row)).lower() for row in rows.tolist())
    assert minimum**2 > p0
    # 在此范围，弱核相对共同参考值的逐层误差上界随截止增大而减小。
    weak_column = lib("full-neutral-free-cutoff")["weak_column"]
    error = (weak_column(minimum**2, p0, p1) - reference).abs_upper()
    weak = reference + arb(0, (error * sum(abs(w).upper() for w in weights)).upper())

    assert arb("1.133895") < newton < arb("1.133905")
    assert arb("-.0204942") < vacuum < arb("-.0204938")
    assert arb("-.07425657719468") < weak < arb("-.07425657719465")
    print("第十三壳盒上 K、W、D 及最小归一化截止下界",
          *[v.str(18) for v in (newton, vacuum, weak, minimum)], flush=True)
    return newton, vacuum, weak


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    fifth, columns, levels = src["verify"](root, return_operator="spectral")
    witness = run_path(str(root / "ch39-six-unrestricted-eighteen.py"))
    readout(root, levels, witness["experiment"](root, fifth, columns, levels))
