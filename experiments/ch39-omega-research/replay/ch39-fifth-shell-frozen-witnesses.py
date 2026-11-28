"""在固定第五壳的四个可行重谱之间，冻结并交换四套正交汇见证。"""
from itertools import product
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
NUMBERS = (2_345_670, 2_345_676, 2_345_684, 2_345_698)


def data(name):
    return json.loads((ROOT/name).read_text())


def witness_moments(load, shell, edges, kernel, layout):
    scale = arb(7)/1000
    def point(i, k):
        a,b,w = edges[i]
        return w*kernel["basis"](scale*((1-k)*a+k*b))
    labels = [arb(k)/11000 for row in layout["base"] for k in row]
    extra = [(None,i,arb(k)/11000) for i,k in layout["extra"]]
    shares, _ = load("one-background-specific-line")["solve_layout"](
        [labels], extra, 0, ["11.5"], [shell], [layout["heavy"]], False, reference_gev=10**12)
    sites = [point(i,k) for i,k in enumerate(labels)]
    measure = sum(sites[1:], sites[0])
    for share, (_,i,k) in zip(shares, extra, strict=True):
        measure += share*(point(i,k)-sites[i])
    return measure


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    exchange = load("fourth-fifth-shell-exchange")
    shell = exchange["FIFTH_SHELL"]
    layouts = [row for row in data("ch39-fourth-fifth-shell-layouts.json") if row["shell"] == shell]
    layouts.extend(data("ch39-fifth-shell-holdout-layouts.json"))
    layouts.append(data("ch39-fifth-shell-quadratic-layout.json"))
    assert tuple(row["heavy"] for row in layouts) == NUMBERS
    ctx = load("shared-large-scale-gap")["load"]()
    fields, beta, kernel = ctx[1], ctx[2], ctx[-2]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    scale = arb(7)/1000
    _, _, bgs = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, NUMBERS, tight_neutral=True)
    assert all(max(z.upper() for _, _, z in bg) < min(arb(1), scale*min(a.lower() for a, _, _ in edges)) for bg in bgs)
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    observables = ((1,0,0,0,0),(0,1,0,0,0),(0,0,1-t,t,-1))
    series = [arb_mat(observables)*kernel["series"](bg, fields) for bg in bgs]
    factor = 2*sum(abs(w) for _, _, w in edges)/(4*arb.pi()*factorial(kernel["ORDER"]+1))
    tails = [[factor*sum(abs(2*row[0]*a)+abs(row[1]*g)+sum(abs(d*q) for d,q in zip(row[2:],bs,strict=True))
                         for (g,a,_),bs in zip(bg,fields,strict=True)) for row in observables] for bg in bgs]
    moments = [witness_moments(load, shell, edges, kernel, layout) for layout in layouts]
    outputs = {}
    for i,j in product(range(len(NUMBERS)), repeat=2):
        if i==j:
            outputs[i,j] = [arb(0)]*3
            continue
        value = (series[j]-series[i])*moments[i]
        outputs[i,j] = [value[k,0]+arb(0,(tails[i][k]+tails[j][k]).upper()) for k in range(3)]
        for delta, bound in zip(outputs[i,j], ("1.99e-8","2.74e-9",".00723"), strict=True):
            assert abs(delta) > arb(bound)
        print(NUMBERS[i], NUMBERS[j], *[x.str(9) for x in outputs[i,j]], flush=True)
    return outputs


if __name__ == "__main__":
    verify()
