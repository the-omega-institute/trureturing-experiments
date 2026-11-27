"""读取按第七、八壳校准的每边一次站位算子在第九壳的真实引力矩。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    ordered = lib("two-shell-ordered-bond-gates")
    box, _ = ordered["experiment"](root, fifth, columns, levels)
    _, _, _, basis, _, weights, heat, field = lib("fifth-seventh-shortest-site-gates")["geometry"](
        root, fifth, columns, levels, number=2_306_888_100_094
    )
    left = basis.conjugate().transpose()
    gravity = lib("shortest-common-complex-gate")["gravity_system"]

    def read(point, values):
        frame, *directions = ordered["transport"](fifth, point)
        projection = (left * frame) * columns
        derivatives = [(left * direction) * columns for direction in directions]
        return gravity(projection, derivatives, levels, weights, heat, field, values)

    mid = [q.mid() for q in box]
    residual, _ = read(mid, True)
    _, jac = read(box, False)
    result = []
    for i, target in enumerate((arb(1), arb("1e-122"))):
        value = target + residual[i, 0] + sum(jac[i, j]*(q - mid[j]) for j, q in enumerate(box))
        print(("牛顿项", "真空项")[i], value.lower().str(16), value.upper().str(16), flush=True)
        result.append(value)
    k, w = result
    assert arb(".122283") < k < arb(".123940")
    assert arb("-.018590") < w < arb("-.018472")
    return result


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
