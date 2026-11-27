"""截止 52 四轮电路的两族五式校准及整盒第六式的严格排除。"""

import argparse
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def certificate(root):
    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    source = load("fifth-full-positive-joint")
    fifth, columns, levels = source["verify"](root, return_operator="spectral")
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    transport = load("two-shell-ordered-bond-gates")["transport"]
    gravity = load("shortest-common-complex-gate")["gravity_system"]
    calibration = load("three-shell-ordered-five")
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_095):
        _, _, _, basis, _, weights, heat, field = geometry(
            root, fifth, columns, levels, number=number
        )
        shells.append((basis.conjugate().transpose(), weights, heat, field))

    def evaluate(angles, directions, *, values=False):
        moved = transport(
            fifth, [angles[j] for j in directions], parameters=angles, active=directions
        )
        residuals, jacobian = [], []
        for left, weights, heat, field in shells:
            projected, *derivatives = [(left * matrix) * columns for matrix in moved]
            response, derivative = gravity(
                projected, derivatives, levels, weights, heat, field, values
            )
            residuals.extend(response.entries())
            jacobian.extend(derivative.tolist())
        return residuals, jacobian

    def verify_family(parameters, active, rows, radius, subdivisions, ratio_limit, bounds):
        center = list(map(arb, parameters))
        free = [i for i in range(12) if i not in active]
        missing, = set(range(6)) - set(rows)
        h, k = arb(radius), arb("1e-6")
        residuals, derivatives = evaluate(center, active, values=True)
        active_matrix = arb_mat([derivatives[i] for i in rows])
        pre = active_matrix.mid().inv()
        assert not pre.det().contains(0)
        identity = arb_mat([[int(i == j) for j in range(5)] for i in range(5)])

        # 固定五个活动角在中心，让另外七个角在其整盒内变化。
        parameter_box = center.copy()
        for j in free:
            parameter_box[j] += arb(0, k)
        _, derivatives = evaluate(parameter_box, free)
        parameter_slopes = (pre * arb_mat([derivatives[i] for i in rows])).tolist()
        offset = (pre * arb_mat([[residuals[i]] for i in rows])).entries()
        center_cost = [
            (v.abs_upper() + k * sum(x.abs_upper() for x in row)) / h
            for v, row in zip(offset, parameter_slopes, strict=True)
        ]

        # 每一子盒均覆盖全部七个自由角，只细分五个活动角。
        largest = arb(0)
        for cell in product(range(subdivisions), repeat=5):
            box = parameter_box.copy()
            for j, index in zip(active, cell, strict=True):
                shift = h * (2 * index + 1 - subdivisions) / subdivisions
                box[j] = center[j] + shift + arb(0, (h / subdivisions).upper())
            _, derivatives = evaluate(box, active)
            defect = identity - pre * arb_mat([derivatives[i] for i in rows])
            ratio = max(
                cost + sum(x.abs_upper() for x in row)
                for cost, row in zip(center_cost, defect.tolist(), strict=True)
            )
            largest = max(largest, ratio)
        assert largest < arb(ratio_limit)

        # 未选读数在整个十二维盒上取均值界，包含每个校准根。
        widths = [h if j in active else k for j in range(12)]
        full_box = [x + arb(0, r) for x, r in zip(center, widths, strict=True)]
        _, derivatives = evaluate(full_box, range(12))
        error = sum(r * d.abs_upper() for r, d in zip(widths, derivatives[missing], strict=True))
        omitted = residuals[missing] + arb(0, error.abs_upper())
        assert arb(bounds[0]) < omitted < arb(bounds[1])
        print("余项、全盒自映比及余项端点", missing, largest.str(12),
              omitted.lower().str(16), omitted.upper().str(16), flush=True)

    verify_family(calibration["ANGLES"], calibration["ACTIVE"], (0, 1, 2, 3, 4),
                  ".00006", 3, ".943", ("-.010113", "-.010064"))
    verify_family(calibration["VACUUM_ANGLES"], calibration["VACUUM_ACTIVE"], (0, 1, 2, 3, 5),
                  ".00005", 2, ".925", ("-.203912", "-.202978"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replay_dir", type=Path, help="钉版 ch39-omega-research/replay 目录")
    certificate(parser.parse_args().replay_dir.resolve())
