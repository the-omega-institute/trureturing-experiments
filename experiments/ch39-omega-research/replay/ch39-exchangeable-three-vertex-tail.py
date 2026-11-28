"""以两个目标的四对小相位扰动检验任意至多三顶点支撑。"""
from itertools import product
from pathlib import Path
from runpy import run_path
from tempfile import TemporaryDirectory
import subprocess
import sys

from flint import acb, arb, arb_mat, ctx


def lower_integer_rows(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    selector = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    ctx.prec = 320
    triple = selector["geometry"](shell, pump)
    labels = list(product(grouped["MASKS"], grouped["MODES"]))
    row = {label: i for i, label in enumerate(labels)}
    geometries = grouped["geometry_classes"]()
    factors = []
    for charge in (-4, 2):
        polys = [selector["squared_factors"](triple(charge, *label)) for label in labels]
        factors.append(arb_mat([[p[k] for k in range(7)] for p in polys]))

    def moments(levels):
        return [sum((e/9-1)**k for e in levels) for k in range(13)]

    frames = list(pump["rotation_bases"](shell))
    adjoint = pump["adjoint"]
    rows = []
    for charge, index in ((-4,2345666), (2,964890)):
        central = moments(light(charge, index))
        base, edges = pump["clock_matrices"](shell, charge)
        phases = shell["geo"]["phases"](acb.zeta_zero(index).imag)
        z = [phases[k]**charge for k in (724,804)]
        for a,b in ((1,0),(-1,0),(0,1),(0,-1)):
            shifts = [acb(0, arb(k)/10**7).exp() for k in (a,b)]
            links = edges[0]*z[0]*shifts[0] + edges[1]*z[1]*shifts[1]
            operator = base + links + adjoint(links)
            energies = [selector["minimum"](adjoint(B)*operator*B) for B in frames]
            delta = [x-y for x,y in zip(central, moments(energies))]
            kernel = arb_mat([[delta[i+j] for j in range(7)] for i in range(7)])
            values = factors[0] * kernel * factors[1].transpose()
            gaps = []
            for pairs in geometries.values():
                for up, down in grouped["MODE_PAIRS"]:
                    gaps.append(sum(values[row[u, up],row[d, down]] for u,d in pairs) / len(pairs))
            rows.append(gaps)
    result = []
    for direction in zip(*rows):
        size = max(abs(z).upper() for z in direction)
        assert size > 0
        result.append([int((z*10**12/size).lower().floor().unique_fmpz()) for z in direction])
    assert len(result) == 1136 and all(len(row) == 8 for row in result)
    return result


def verify(replay):
    matrix = lower_integer_rows(replay)
    source = Path(__file__).with_suffix(".cpp")
    serialized = "".join(" ".join(map(str, row))+"\n" for row in matrix)
    with TemporaryDirectory() as work:
        program = Path(work) / "three-vertex-tail"
        subprocess.run(["c++","-std=c++17","-O3",str(source),"-o",str(program)],check=True)
        subprocess.run([str(program)],input=serialized,text=True,check=True)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent)
