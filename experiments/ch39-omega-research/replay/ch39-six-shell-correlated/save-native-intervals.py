"""Save exact Arb input enclosures for local certificate development."""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def encode(real):
    return [[int(a) for a in point.man_exp()] for point in (real.mid(), real.rad())]


def matrix(a):
    return [[[encode(z.real), encode(z.imag)] for z in row] for row in a.tolist()]


def main(root):
    ctx.prec = 192
    pinned = root / "pinned"
    fifth = run_path(str(pinned / "ch39-fifth-full-positive-joint.py"))
    basis, columns, levels = fifth["verify"](pinned, return_operator="spectral")
    geometry = run_path(str(pinned / "ch39-fifth-seventh-shortest-site-gates.py"))["geometry"]
    native = {
        "basis": matrix(basis),
        "columns": matrix(columns),
        "levels": list(map(encode, levels)),
        "shells": [],
    }
    for n in range(2_306_888_100_092, 2_306_888_100_098):
        _, _, _, shell_basis, _, weights, _, field = geometry(
            pinned, basis, columns, levels, number=n
        )
        left = shell_basis.conjugate().transpose()
        native["shells"].append(
            [matrix(left), list(map(encode, weights)),
             [[encode(arb(g)), encode(a), encode(m)] for g, a, m in field]]
        )
        print(n, flush=True)
    target = root / "native-intervals.json"
    pending = target.with_suffix(".tmp")
    pending.write_text(json.dumps(native))
    pending.replace(target)


if __name__ == "__main__":
    main(Path(__file__).resolve().parent)
