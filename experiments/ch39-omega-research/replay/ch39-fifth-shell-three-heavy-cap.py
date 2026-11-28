"""只用第一个重谱，也能给任意自由分配的正交汇能标定上界。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
DIRECTION = (268787, -330912, -5772331377, 323745708, -7373213999)
LOWER = 5_126_800_000_000_000


def verify():
    four = run_path(str(ROOT/"ch39-four-background-crossing.py"))
    load = four["load"]
    core = load("shared-large-scale-gap")
    ctx = core["load"]()
    fields, beta, kernel = ctx[1], ctx[2], ctx[-2]
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    positive = load("fifth-shell-common-positive")
    _, lib, (bg,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT, positive["HEAVY"][:1], tight_neutral=True)
    assert max(z.upper() for _, _, z in bg) < min(arb(1), min(a.lower() for a, _, _ in edges)*arb(7)/1000)
    single = ([(edges, kernel["series"](bg, fields), bg)], *ctx[1:])
    model = lib["chiral"]["source"]["mean"]["inputs"]()[0]
    colour = arb(DIRECTION[4])/10**10
    assert colour < 0
    vector = [10*d for d in DIRECTION[:4]]
    groups, constant, _ = core["project_direction"](
        single, vector, kernel_ratio=1, low_colours=[colour],
        common_log=(arb(10**16)/model["mz"]).log())
    bound = four["lower_projection"](single, groups, constant, 0)
    slope = -(DIRECTION[2]*(beta[0]-beta[2])+DIRECTION[3]*(beta[1]-beta[2]))/10**10
    assert slope > 0
    upper = 10**16*(-bound/slope).exp()
    print("第一个重谱全部谱线独立时，正交汇能标的严格上界 GeV", upper.str(14), flush=True)
    assert upper < 5_127_109_000_000_000
    layouts = json.loads((ROOT/"ch39-fifth-shell-edge-layouts.json").read_text())
    assert sorted(item["line"] for item in layouts) == [285, 323]
    for layout in layouts:
        positive["solve_witness"](load, shell, layout["line"], layout=layout,
                                  energy=LOWER, minimum=".00003")


if __name__ == "__main__":
    verify()
