"""Certified four-shell targets for the fixed last-round sign and varied shell-ten readings."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


ANGLES = (
    "1.024629488089784", "-0.08154296875", "-2.237196091186837", "2.134765625",
    "-2.794677734375", "3.068603515625", "-0.9077401330013921", "2.9622941968913445",
    "1.0860539997121788", "-2.140460614670698", "3.038818380548489", "-1.1722108480592708",
)
ACTIVE = (0, 2, 6, 7, 8, 9, 10, 11)
RADII = ("0.000050", "0.000076", "0.000080", "0.000083",
         "0.000022", "0.000091", "0.000025", "0.000056")


def verify(root):
    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    source = load("fifth-full-positive-joint")
    fifth, columns, levels = source["verify"](root, return_operator="spectral")
    station = fifth * columns
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    ordered = load("two-shell-ordered-bond-gates")
    gravity = load("shortest-common-complex-gate")["gravity_system"]
    centered = run_path(str(Path(__file__).with_name("ch39-three-shell-final-sign.py")))["centered_residual"]
    layers = ordered["LAYERS"]
    assert [k for k, _ in layers] == [0, 1, 2, 2]
    endpoints = [site for edge in layers[-1][1] for site in edge]
    assert len(set(endpoints)) == len(endpoints) == 6
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_096):
        _, _, _, basis, _, weights, heat, fields = geometry(root, fifth, columns, levels, number=number)
        left = basis.conjugate().transpose()
        # Absorb the last-round signs into the basis before projecting any derivative.
        for site in endpoints:
            for row in range(left.nrows()):
                left[row, site] = -left[row, site]
        shells.append((left, weights, heat, fields))
    center = [arb(ANGLES[k]) for k in ACTIVE]
    radii = list(map(arb, RADII))
    transport = ordered["transport"](station, center, parameters=ANGLES, active=ACTIVE)
    values, jacobian = [], []
    for left, weights, heat, fields in shells:
        moved, *slopes = [left * m for m in transport]
        values.extend(centered(moved, levels, weights, heat, fields).tolist())
        _, jac = gravity(moved, slopes, levels, weights, heat, fields, False)
        jacobian.extend(jac.tolist())
    jacobian = arb_mat(jacobian)
    pre = jacobian.mid().inv().mid()
    assert not pre.det().contains(0)
    offset = pre * arb_mat(values)
    ident = arb_mat([[int(i == j) for j in range(8)] for i in range(8)])
    defect = ident - pre * jacobian

    # This pinned Hessian routine differentiates six axes at once. Each call
    # still evaluates on the full eight-axis box, including the other two axes.
    spec = spec_from_file_location("shell_curvature", root / "ch39-three-shell-gravity-radius56.py")
    second = module_from_spec(spec)
    spec.loader.exec_module(second)
    params = list(map(arb, ANGLES))
    for k, c, r in zip(ACTIVE, center, radii, strict=True):
        params[k] = c + arb(0, r.upper())
    second.PARAMETERS = params
    bounds = {}
    for group in ((0, 1, 2, 3, 4, 5), (0, 1, 2, 3, 6, 7), (2, 3, 4, 5, 6, 7)):
        second.ACTIVE = [ACTIVE[i] for i in group]
        box = [params[k] for k in second.ACTIVE]
        hessian, pairs = second.hessian(station, layers, box, ordered["block"], shells, levels)
        weighted = pre * hessian
        for col, (a, b) in enumerate(pairs):
            bounds[group[a], group[b]] = [weighted[i, col].abs_upper() for i in range(8)]
    assert set(bounds) == {(a, b) for a in range(8) for b in range(a + 1)}

    # Weighted sup norm on the angle box. The Hessian bounds enclose the
    # variation of the first derivative over every segment from the center.
    norms, ratios = [], []
    for i, r in enumerate(radii):
        linear = sum(defect[i, j].abs_upper() * radii[j] for j in range(8))
        variation = sum(row[i] * radii[a] * radii[b] * (1 if a == b else 2)
                        for (a, b), row in bounds.items())
        slope = (linear + variation) / r
        shift = arb(".000004") * pre[i, 6].abs_upper() + arb(".0000005") * pre[i, 7].abs_upper()
        norms.append(slope)
        ratios.append(slope + (offset[i, 0].abs_upper() + shift) / r)
    print("Weighted contraction bound", max(norms).str(14))
    print("Uniform self-map bound", max(ratios).str(14))
    assert all(z < arb(".023") for z in norms)
    assert all(z < 1 for z in ratios)


if __name__ == "__main__":
    ctx.prec = 192
    verify(Path(sys.argv[1]).resolve())
