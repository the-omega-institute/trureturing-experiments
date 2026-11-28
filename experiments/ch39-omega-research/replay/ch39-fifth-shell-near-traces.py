"""四十八个新重谱中，检验全谱迹预报量最小的五对。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
NUMBERS = tuple(range(2_345_772, 2_345_867, 2))


def nearest(load, numbers):
    geometry = load("three-ray-chiral-boundary")
    trace = load("fifth-shell-closed-walk-vacuum")["full_trace"]
    ctx.prec = 320
    values = {n:trace(geometry,n) for n in numbers}
    forecasts = {(a,b):(values[b]-values[a])/(-arb(103)/3) for a,b in combinations(numbers,2)}
    ranking = sorted(forecasts, key=lambda p: abs(forecasts[p]).mid())
    assert len(ranking)==1128
    for i,p in enumerate(ranking[:5]):
        assert all(abs(forecasts[p]) < abs(forecasts[q]) for q in ranking[i+1:])
    return {p:forecasts[p] for p in ranking[:5]}


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    forecasts = nearest(load,NUMBERS)
    pairs = list(forecasts)
    numbers = sorted({n for p in pairs for n in p})
    intervals = load("fifth-shell-vacuum-order")["vacuum_intervals"]
    wrong = []
    for a,b,lead,change,_ in intervals(load,numbers,selected_pairs=pairs):
        pred = forecasts[a,b]
        assert not (pred*change).contains(0) and lead*change>0
        if pred*change<0:
            wrong.append((a,b))
        print(a,b,"全谱预报、重谱主项、真空差下上端",pred.str(11),lead.str(11),
              change.lower().str(11),change.upper().str(11),flush=True)

    assert wrong == [(2_345_776,2_345_794),(2_345_778,2_345_852)]

    # 两种明确分配，逐条让正负权分别取相反端点，再以相同份额混合。
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT,load("native-rh-ray-pvm")["shell_weights"],shell)
    sources = sorted({a for a,b in pairs})
    _,lib,bgs = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,sources,tight_neutral=True)
    heat = lib["expint"]
    scale = arb(7)/1000
    for source,rows in zip(sources,bgs,strict=True):
        low,high = arb(0),arb(0)
        for a,b,w in edges:
            def newton(x):
                return x*sum(alpha*heat(2,z/x) for _,alpha,z in rows)/(2*arb.pi())
            left,right = w*newton(scale*a),w*newton(scale*b)
            low += left.min(right)
            high += left.max(right)
        assert low < 1 < high
        print(source,"两种端点分配的源牛顿矩",low.str(9),high.str(9),flush=True)


if __name__=="__main__":
    verify()
