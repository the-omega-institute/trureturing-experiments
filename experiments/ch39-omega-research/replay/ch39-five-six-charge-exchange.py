"""第五壳固定权、同一重谱下，独立交换荷 −4 与荷 +2 的归一化本征值。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
SHELL = {5:2306888100089,6:2306888100090}
HEAVY = 2345676


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def primal():
    layouts = json.loads((ROOT/"ch39-five-six-charge-exchange-layouts.json").read_text())
    assert [item["sectors"] for item in layouts]==[[5,5],[5,6],[6,5]]
    old = json.loads((ROOT/"ch39-fifth-shell-holdout-layouts.json").read_text())
    native, = [item for item in old if item["heavy"]==HEAVY]
    solve = load("one-background-specific-line")["solve_layout"]
    cases = [(native,5,5,7000,"11.71")]
    for layout,minimum in zip(layouts,("11.77","11.74","11.75"),strict=True):
        cases.append((layout,*layout["sectors"],5000,minimum))
    for layout,minus,plus,scale,minimum in cases:
        labels = [arb(k)/11000 for row in layout["base"] for k in row]
        extra = [(None,i,arb(k)/11000) for i,k in layout["extra"]]
        print("荷 −4 壳、荷 +2 壳、标尺（百万分之一）",minus,plus,scale,flush=True)
        solve([labels],extra,set(),[minimum],[SHELL[minus]],[HEAVY],False,reference_gev=10**12,
              scale_millionths=[scale],weight_codes=[SHELL[5]],positive_codes=[SHELL[plus]])


def dual():
    _,rows,_,model,initial,beta,fields,_,_ = load("sixth-all-scales")["inputs"]()
    kernel = load("continuous-split-roots")
    series = kernel["series"](rows,fields)
    paired = load("shared-twelve-node-gauge-gap")["paired_lines"]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    project = load("three-background-calibration")["projection"]
    cases = [(5,6,7000,(-59118461522,-100000000000,34475326043,-100000000000,0)),
             (6,5,7000,(-7939401030,-100000000000,2902541510,-8418081017,0)),
             (6,6,7000,(-57642435613,-100000000000,34475176484,-100000000000,0)),
             (6,6,5000,(-5166543295,-100000000000,541403203,-1568992051,0))]
    for (minus,plus,scale,direction),lower in zip(cases,(".47179",".78681","1.3769",".42344"),strict=True):
        edges = paired(ROOT,shell,SHELL[minus],positive_code=SHELL[plus],weight_code=SHELL[5])
        bound,_ = project([(edges,series,rows)],fields,beta,initial,model,kernel,[direction],None,set(),
                          scale_millionths=[scale])
        norm = sum(abs(arb(v)/10**11) for v in direction[2:4])
        print("荷 −4 壳、荷 +2 壳、标尺与规范差下界",minus,plus,scale,(bound/norm).str(13),flush=True)
        assert bound/norm>arb(lower)


if __name__=="__main__":
    primal()
    dual()
