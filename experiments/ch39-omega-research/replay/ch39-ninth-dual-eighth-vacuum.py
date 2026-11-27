"""同一组四轮十二角同时校准第八壳真空矩与第九壳两个原生引力矩。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


ANGLES = (".359375", "1.390625", "3", "1.5625", "-1.16482345221599", "-7.359375",
          "7.375", "-8.25", "1.07247168750305", "-4.5", "-5.453125", "-.799293898988015")
ACTIVE = (4, 8, 11)
RADIUS = ".00018"


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    transport = lib("two-shell-ordered-bond-gates")["transport"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    station = fifth * columns
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_095):
        _, _, _, basis, _, weights, heat, field = geometry(root, fifth, columns, levels, number=number)
        shells.append((basis.conjugate().transpose(), weights, heat, field))

    def system(point, data, values):
        transported = transport(station, point, parameters=ANGLES, active=ACTIVE)
        residuals, jacobian = [], []
        for left, weights, heat, field in data:
            projection, *directions = [left * m for m in transported]
            f, j = gravity(projection, directions, levels, weights, heat, field, values)
            residuals.extend(f.entries())
            jacobian.extend(j.tolist())
        return residuals, jacobian

    center = [arb(ANGLES[i]) for i in ACTIVE]
    radius = arb(RADIUS)
    # 第八、九壳四行依次为 K8、W8、K9、W9；后三行参与校准。
    f, j = system(center, shells[1:], True)
    pre = arb_mat(j[1:]).mid().inv()
    assert not pre.det().contains(0)
    eye = arb_mat([[int(i == j) for j in ACTIVE] for i in ACTIVE])
    offset = max(v.abs_upper() for v in (pre * arb_mat([[v] for v in f[1:]])).entries())
    defect = arb(0)
    for cell in product(range(4), repeat=3):
        sub = [c + radius*(2*k-3)/4 + arb(0, (radius/4).upper()) for c, k in zip(center, cell, strict=True)]
        _, directions = system(sub, shells[1:], False)
        error = eye - pre * arb_mat(directions[1:])
        defect = max(defect, max(sum(v.abs_upper() for v in row) for row in error.tolist()))
    ratio = offset/radius + defect
    gain = max(sum(v.abs_upper() for v in row) for row in pre.tolist())
    tolerance = radius*(1-ratio)/gain
    print("中心位移、导数偏差、自映比、三目标可独立扰动的共同半径下界",
          *[v.str(12) for v in (offset, defect, ratio, tolerance)], flush=True)
    assert ratio < arb(".913") and tolerance > arb("2.01e-7")

    box = [c + arb(0, radius.upper()) for c in center]
    f, _ = system(center, shells, True)
    _, j = system(box, shells, False)
    # 公共根处用整盒导数的均值界读出未被校准的三项。
    values = [f[i] + (arb(1) if i % 2 == 0 else arb("1e-122"))
              + sum(g*arb(0, radius.upper()) for g in j[i]) for i in range(3)]
    print("第七壳牛顿矩、第七壳真空矩、第八壳牛顿矩", *[v.str(12) for v in values], flush=True)
    for value, lo, hi in zip(values, ("1.327", ".058", "1.548"), ("1.333", ".060", "1.561"), strict=True):
        assert arb(lo) < value < arb(hi)
    return box, ratio, tolerance, values


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
