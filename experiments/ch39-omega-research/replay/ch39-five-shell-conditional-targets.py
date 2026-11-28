"""Joint shell-eleven target rectangles conditional on the shell-seven-to-ten calibration."""

from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def verify(root):
    here = Path(__file__).resolve().parent

    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    fifth, columns, levels = load("fifth-full-positive-joint")["verify"](
        root, return_operator="spectral"
    )
    station = fifth * columns
    ordered = load("two-shell-ordered-bond-gates")
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = load("shortest-common-complex-gate")["gravity_system"]
    residual = run_path(str(here / "ch39-three-shell-final-sign.py"))[
        "centered_residual"
    ]
    spec = spec_from_file_location(
        "shell_second", root / "ch39-three-shell-gravity-radius56.py"
    )
    second = module_from_spec(spec)
    spec.loader.exec_module(second)
    layers = ordered["LAYERS"]
    assert [family for family, _ in layers] == [0, 1, 2, 2]
    endpoints = [s for edge in layers[-1][1] for s in edge]
    assert len(endpoints) == len(set(endpoints)) == 6
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_097):
        _, _, _, basis, _, weights, heat, field = geometry(
            root, fifth, columns, levels, number=number
        )
        left = basis.conjugate().transpose()
        for site in endpoints:
            for row in range(left.nrows()):
                left[row, site] = -left[row, site]
        shells.append((left, weights, heat, field))

    charts = json.loads(Path(__file__).with_suffix(".json").read_text())
    assert sorted(tuple(c["target"]) for c in charts) == [
        ("0.4", "0.05"),
        ("0.4", "0.06"),
        ("0.44", "0.05"),
        ("0.44", "0.06"),
    ]
    eye = arb_mat([[int(i == j) for j in range(10)] for i in range(10)])
    for chart in charts:
        active = chart["active"]
        params = list(map(arb, chart["angles"]))
        radii = list(map(arb, chart["radii"]))
        assert len(params) == 12 and len(set(active)) == len(active) == len(radii) == 10
        assert min(radii) > 0
        k, w = map(arb, chart["target"])
        center = [params[i] for i in active]
        values, jacobian = [], []
        moved = ordered["transport"](station, center, parameters=params, active=active)
        for left, weights, heat, field in shells:
            projected, *slopes = [left * m for m in moved]
            values.extend(residual(projected, levels, weights, heat, field).tolist())
            _, jac = gravity(projected, slopes, levels, weights, heat, field, False)
            jacobian.extend(jac.tolist())
        # The reused residual subtracts the original two targets at every shell.
        values[8][0] += 1 - k
        values[9][0] += arb("1e-122") - w
        jacobian = arb_mat(jacobian)
        pre = jacobian.mid().inv().mid()
        assert not pre.det().contains(0)
        offset = pre * arb_mat(values)
        defect = eye - pre * jacobian
        box = [c + arb(0, r.upper()) for c, r in zip(center, radii, strict=True)]
        for i, z in zip(active, box, strict=True):
            params[i] = z
        second.PARAMETERS = params
        curvature = {}
        # The pinned routine differentiates six axes while bounding all ten over the box.
        for group in (
            (0, 1, 2, 3, 4, 5),
            (4, 5, 6, 7, 8, 9),
            (0, 1, 6, 7, 8, 9),
            (2, 3, 6, 7, 8, 9),
        ):
            second.ACTIVE = [active[i] for i in group]
            hessian, pairs = second.hessian(
                station,
                layers,
                [box[i] for i in group],
                ordered["block"],
                shells,
                levels,
            )
            bounded = pre * hessian
            for col, (a, b) in enumerate(pairs):
                curvature[group[a], group[b]] = [
                    bounded[i, col].abs_upper() for i in range(10)
                ]
        assert set(curvature) == {(a, b) for a in range(10) for b in range(a + 1)}
        ratios = []
        for i, r in enumerate(radii):
            linear = sum(defect[i, j].abs_upper() * s for j, s in enumerate(radii))
            curve = sum(
                row[i] * radii[a] * radii[b] * (1 if a == b else 2)
                for (a, b), row in curvature.items()
            )
            shift = pre[i, 8].abs_upper() * arb(".000003") + pre[
                i, 9
            ].abs_upper() * arb(".000000375")
            ratios.append((offset[i, 0].abs_upper() + shift + linear + curve) / r)
        assert all(z < 1 for z in ratios), (chart["target"], ratios)
        print(
            "Joint targets and uniform self-map bound",
            *chart["target"],
            max(ratios).str(10),
            flush=True,
        )


if __name__ == "__main__":
    ctx.prec = 192
    verify(Path(sys.argv[1]).resolve())
