"""Five native Newton/vacuum targets with one exceptional bond and a shell-twelve holdout."""

import json
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def held_out_readings(model, shell, layers, chart, interval_angles):
    station, levels, ordered, gravity, residual, _, _ = model
    left, weights, heat, field = shell
    active = chart["active"]
    params = list(map(arb, chart["angles"]))
    radii = list(map(arb, chart["radii"]))
    values = ordered["transport"](station, [], parameters=params, active=[], layers=layers)
    offset = residual(left * values[0], levels, weights, heat, field)
    box = [interval_angles[i] for i in active]
    jets = ordered["transport"](station, box, parameters=params, active=active, layers=layers)
    projected, *derivatives = [left * m for m in jets]
    _, slopes = gravity(projected, derivatives, levels, weights, heat, field, False)
    bounds = []
    for i, target in enumerate((arb(1), arb("1e-122"))):
        error = sum(slopes[i, j].abs_upper() * radius for j, radius in enumerate(radii))
        bounds.append(target + offset[i, 0] + arb(0, error.upper()))
    return bounds


def verify(root):
    shared = run_path(str(Path(__file__).with_name("ch39-five-shell-independent-fourth.py")))
    model = shared["native_model"](root)
    _, _, ordered, _, _, _, shell_at = model
    shells = [shell_at(n) for n in range(2_306_888_100_092, 2_306_888_100_097)]
    rounds = ordered["LAYERS"]
    assert [g for g, _ in rounds] == [0, 1, 2, 2]
    edges = [(i, edge) for i, (_, pairs) in enumerate(rounds) for edge in pairs]
    charts = json.loads(Path(__file__).with_suffix(".json").read_text())
    assert [c["split"] for c in charts] == [9, 4, 6, 8, 10]
    assert [edges[c["split"]][1] for c in charts] == [(6, 13), (1, 3), (6, 9), (8, 16), (7, 10)]
    for chart in charts:
        assert chart["newton"] == "1"
        split = chart["split"]
        layers = [(4 if i == split else family, (edge,)) for i, (family, edge) in enumerate(edges)]
        angles, ratio = shared["certify_chart"](model, shells, layers, chart)
        print("Independent edge and maximum contraction/self-map ratio",
              edges[split][1], ratio.str(10), flush=True)
        if split == 9:
            bounds = held_out_readings(model, shell_at(2_306_888_100_097), layers, chart, angles)
            for reading, lower, upper in zip(bounds, ("0.170", "0.1030"), ("0.172", "0.1033"), strict=True):
                assert arb(lower) < reading < arb(upper), reading
            print("Held-out Newton and vacuum intervals", *[b.str(10) for b in bounds], flush=True)


if __name__ == "__main__":
    ctx.prec = 192
    verify(Path(sys.argv[1]).resolve())
