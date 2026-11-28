"""固定第五壳和三个重谱，分类 10¹² GeV 下全部 361 条可独立谱线。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
MEETING = 10**12


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    positive = load("fifth-shell-common-positive")
    kernel = load("continuous-split-roots")
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    _, lib, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT, positive["HEAVY"], tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0, 0, 0), *load("signed-gravity-unification")["gauge_weights"](model)]
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    sectors = [(edges, kernel["series"](bg, fields), bg) for bg in spectra]
    build = load("three-background-calibration")["projection_terms"]
    minimum = load("shared-continuous-position-gap")["projected_minimum"]
    directions = json.loads((ROOT/"ch39-fifth-shell-line-directions.json").read_text())
    remaining = set(range(361))
    movable = remaining.copy()
    for flat in directions:
        # 对每个背景依次取 K-1、W-目标、第一与第三支之差、第二与第三支之差。
        vectors = [[*flat[j:j+4], 0] for j in (0, 4, 8)]
        base, slope, groups = build(sectors, fields, beta, initial, model, kernel, vectors, MEETING, None)
        norm = sum(abs(arb(k)/10**11) for row in vectors for k in row[2:4])
        def lower(terms):
            return minimum(kernel["basis"], kernel["S"], terms, arb("1e-11"))
        shared = []
        for row in groups:
            # 三处是完全相同的壳和标尺，分布共用时先合并热核系数。
            left, span, coeff = row[0]
            coeff = coeff + row[1][2] + row[2][2]
            shared.append(lower([(left, span, coeff)]))
        total = base + sum(shared)
        bounds = {i: (total-shared[i]+sum(lower([term]) for term in groups[i]))/norm
                  for i in remaining | movable}
        eliminated = {i for i in remaining if bounds[i] > arb(".00105")}
        assert eliminated
        print("本步排除条数、最小规范差下界", len(eliminated),
              min(bounds[i].lower() for i in eliminated).str(12), flush=True)
        remaining -= eliminated
        shift = min((slope*(energy/MEETING).log()).lower() for energy in (model["mz"], arb("1.22089e19")))
        wide = {i: bounds[i]+shift/norm for i in movable}
        excluded = {i for i in movable if wide[i] > arb(".00022")}
        print("共同能标全域本步排除条数", len(excluded), flush=True)
        movable -= excluded
    assert remaining == movable == {285, 323}
    for line in sorted(remaining):
        positive["solve_witness"](load, shell, line)
    layouts = json.loads((ROOT/"ch39-fifth-shell-common-energy-layouts.json").read_text())
    assert {(entry["line"], entry["energy"]) for entry in layouts} == {
        (line, energy) for line in movable for energy in (92, 5*10**15)}
    for entry in layouts:
        print("交汇能标、独立线", entry["energy"], entry["line"], flush=True)
        positive["solve_witness"](load, shell, entry["line"], layout=entry,
                                  energy=entry["energy"], minimum=".033")


if __name__ == "__main__":
    verify()
