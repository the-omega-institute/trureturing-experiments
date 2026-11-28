"""Five exact vacuum targets and four Newton targets with an independent fourth gate."""

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
    original = ordered["LAYERS"]
    assert [family for family, _ in original] == [0, 1, 2, 2]
    layers = (*original[:3], (3, original[3][1]))
    # Sites 9 and 13 are second endpoints exclusive to rounds three and four.
    for leaf, expected in ((9, 2), (13, 3)):
        hits = [(j, pair.index(leaf)) for j, (_, edges) in enumerate(layers)
                for pair in edges if leaf in pair]
        assert hits == [(expected, 1)]
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_097):
        _, _, _, basis, _, weights, heat, field = geometry(
            root, fifth, columns, levels, number=number
        )
        left = basis.conjugate().transpose()
        shells.append((left, weights, heat, field))

    charts = json.loads(Path(__file__).with_suffix(".json").read_text())
    assert [c["newton"] for c in charts] == ["0.7", "0.72"]
    eye = arb_mat([[int(i == j) for j in range(10)] for i in range(10)])
    for chart in charts:
        active = chart["active"]
        params = list(map(arb, chart["angles"]))
        radii = list(map(arb, chart["radii"]))
        assert len(params) == 16 and len(set(active)) == len(active) == len(radii) == 10
        assert all(0 <= j < 16 for j in active)
        assert min(radii) > 0
        target = arb(chart["newton"])
        center = [params[i] for i in active]
        values, jacobian = [], []
        moved = ordered["transport"](
            station, center, parameters=params, active=active, layers=layers
        )
        for left, weights, heat, field in shells:
            projected, *slopes = [left * m for m in moved]
            values.extend(residual(projected, levels, weights, heat, field).tolist())
            _, jac = gravity(projected, slopes, levels, weights, heat, field, False)
            jacobian.extend(jac.tolist())
        # The reused residual subtracts the original two targets at every shell.
        values[8][0] += 1 - target
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
            # Taylor's self-map bound has a half; the derivative bound uses the full Hessian.
            self_map = offset[i, 0].abs_upper() + linear + curve / 2
            ratios.append(max(self_map, linear + curve) / r)
        assert all(z < 1 for z in ratios), (target, ratios)
        distance = abs(abs(params[8].cos()) - abs(params[12].cos())) / 2
        assert distance > arb("0.01")
        print("Newton target, maximum ratio, phase-only distance floor",
              target.str(4), max(ratios).str(10), distance.lower().str(8), flush=True)



if __name__ == "__main__":
    ctx.prec = 192
    verify(Path(sys.argv[1]).resolve())
