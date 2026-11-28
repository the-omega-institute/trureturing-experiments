"""Bound input and angle errors separately for the six-shell calibration."""
import json
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def real(pair):
    return arb(tuple(pair[0]), tuple(pair[1]))


def matrix(items):
    return acb_mat([[acb(real(z[0]), real(z[1])) for z in row] for row in items])


def verify(root):
    ctx.prec = 256
    pinned = root / "pinned"
    native = json.loads((root / "native-intervals.json").read_text())
    candidate = json.loads((root / "candidate.json").read_text())
    basis, columns = [matrix(native[key]) for key in ("basis", "columns")]
    levels = list(map(real, native["levels"]))
    ordered = run_path(str(pinned / "ch39-two-shell-ordered-bond-gates.py"))
    gravity = run_path(str(pinned / "ch39-shortest-common-complex-gate.py"))["gravity_system"]
    residual = run_path(str(root / "ch39-three-shell-final-sign.py"))["centered_residual"]
    heat = run_path(str(pinned / "ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    second = run_path(str(pinned / "ch39-three-shell-gravity-radius56.py"))["hessian"]
    shell_inputs = [(matrix(left), list(map(real, weights)), heat,
                     [tuple(map(real, triple)) for triple in field])
                    for left, weights, field in native["shells"]]

    original = [(i, edge) for i, (_, pairs) in enumerate(ordered["LAYERS"]) for edge in pairs]
    split, extra = candidate["split"], candidate["extra"]
    layers = [(4 if j==split else 5 if j==extra else family, (edge,))
              for j, (family, edge) in enumerate(original)]
    parameters = list(map(arb, candidate["angles"]))
    active = candidate["active"]
    center = [parameters[i] for i in active]
    mid_columns = acb_mat([[acb(z.real.mid(),z.imag.mid()) for z in row] for row in columns.tolist()])
    station = basis * columns
    station_mid = basis * mid_columns
    eye = arb_mat([[int(i == j) for j in range(12)] for i in range(12)])

    def evaluate(source, point, want_values):
        entries = ordered["transport"](source, point, parameters=parameters, active=active, layers=layers)
        values, derivatives = [], []
        for left, weights, heat, field in shell_inputs:
            p, *slopes = [left*m for m in entries]
            if want_values:
                values.extend(residual(p, levels, weights, heat, field).tolist())
            _, partial = gravity(p, slopes, levels, weights, heat, field, False)
            derivatives.extend(partial.tolist())
        return arb_mat(values) if want_values else None, arb_mat(derivatives)

    for _ in range(3):
        value, jac = evaluate(station_mid, center, True)
        correction = jac.mid().solve(value.mid())
        center = [(x-correction[j,0]).mid() for j,x in enumerate(center)]
    centers = [x.str(40,radius=False) for x in center]
    center = list(map(arb, centers))
    value, jac = evaluate(station_mid, center, True)
    pre = jac.mid().inv().mid()
    assert not pre.det().contains(0)
    source_radii = [
        z.real.rad() if part == 0 else z.imag.rad()
        for z in columns.entries() for part in range(2)
    ]
    flat_parameters = parameters.copy()
    for index, x in zip(active, center, strict=True):
        flat_parameters[index] = x
    moved_basis, = ordered["transport"](
        basis, [], parameters=flat_parameters, active=[], layers=layers
    )
    source_derivatives = []
    for left, weights, heat, field in shell_inputs:
        transform = left*moved_basis
        project = transform*columns
        # The shared verifier bounds exposures on the physical amplitude range.
        # Check it also contains the entire midpoint-to-source interpolation.
        low, heavy, high = levels
        for row in project.tolist():
            amplitude = low + (heavy-low)*abs(row[0])**2 + (high-low)*abs(row[1])**2
            assert amplitude.upper() < high.upper()
        directions = []
        for i in range(19):
            for j in range(2):
                for scale in (acb(1), acb(0,1)):
                    direction = acb_mat(19,2)
                    for k in range(19):
                        direction[k,j] = scale*transform[k,i]
                    directions.append(direction)
        _, src_jac = gravity(project, directions, levels, weights, heat, field, False)
        source_derivatives.extend(src_jac.tolist())
    pre_source = (pre * arb_mat(source_derivatives)).tolist()
    mid_offset = [x.abs_upper() for x in (pre*value).entries()]
    offset = [a+sum(d.abs_upper()*r for d,r in zip(row,source_radii,strict=True))
              for a,row in zip(mid_offset,pre_source,strict=True)]
    print("Midpoint and total preconditioned offsets",
          max(mid_offset).str(10), max(offset).str(10), flush=True)
    _, jacobian = evaluate(station, center, False)
    defects = (eye-pre*jacobian).tolist()
    for factor in ("1.1", "1.4", "2", "4", "1.02"):
        radii = [(arb(factor)*max(d,arb("1e-12"))).upper() for d in offset]
        box = [x+arb(0,r) for x,r in zip(center,radii,strict=True)]
        parameters_box = parameters.copy()
        for j,x in zip(active,box,strict=True):
            parameters_box[j] = x
        second.__globals__["PARAMETERS"] = parameters_box
        curvatures = {}
        for a,b in combinations(range(4),2):
            group = [i for cluster in (a,b) for i in range(3*cluster,3*cluster+3)]
            second.__globals__["ACTIVE"] = [active[j] for j in group]
            hf, pairs = second(station, layers, [box[j] for j in group], ordered["block"],
                               shell_inputs, levels)
            pre_hessian = pre * hf
            for col,(i,j) in enumerate(pairs):
                curvatures[group[i],group[j]] = [pre_hessian[k,col].abs_upper() for k in range(12)]
        assert len(curvatures)==78
        linear = [sum(a.abs_upper()*d for a,d in zip(row,radii,strict=True)) for row in defects]
        quadratic = [sum(terms[k]*radii[i]*radii[j]*(1 if i==j else 2)
                         for (i,j),terms in curvatures.items()) for k in range(12)]
        derivative = max((a+b)/r for a,b,r in zip(linear,quadratic,radii,strict=True))
        self_map = max((x+a+b/2)/r for x,a,b,r in zip(offset,linear,quadratic,radii,strict=True))
        print("Half-width / derivative / self-map", max(radii).str(10), derivative.str(10),
              self_map.str(10), flush=True)
        if max(derivative,self_map) < 1:
            result = {**candidate,"active_centers":centers,
                      "radii": [r.str(40, radius=False) for r in radii]}
            (root / "improved-certificate.json").write_text(json.dumps(result,indent=2)+"\n")
            print("Six-shell common root certified.", flush=True)
            return
    raise ArithmeticError("Interval refinement did not produce a contraction.")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
