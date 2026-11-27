"""按预登记协议在第十壳读取十八条有向边的三壳共同算子。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    tri = lib("three-shell-gravity-radius56")
    ordered = lib("two-shell-ordered-bond-gates")
    box = tri["experiment"](root, fifth, columns, levels)
    layers = ordered["LAYERS"] + tri["EXTRA_LAYERS"]
    _, _, _, basis, _, weights, heat, field = lib("fifth-seventh-shortest-site-gates")["geometry"](
        root, fifth, columns, levels, number=2_306_888_100_095
    )
    gravity = lib("shortest-common-complex-gate")["gravity_system"]

    def quantities(point, values):
        out = ordered["transport"](fifth*columns, point, parameters=tri["PARAMETERS"],
                                   active=tri["ACTIVE"], layers=layers)
        projection, *directions = [basis.conjugate().transpose()*m for m in out]
        return gravity(projection, directions, levels, weights, heat, field, values)

    center = [q.mid() for q in box]
    f, _ = quantities(center, True)
    _, j = quantities(box, False)
    readings = [target+f[i,0]+sum(j[i,k]*(q-center[k]) for k,q in enumerate(box))
                for i,target in enumerate((arb(1),arb("1e-122")))]
    for name, value in zip(("牛顿项", "真空项"), readings, strict=True):
        print("第十壳", name, value.lower().str(15), value.upper().str(15), flush=True)
    k, w = readings
    assert arb("1.07020") < k < arb("1.07540")
    assert arb("-.002820") < w < arb("-.002444")
    return readings


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
