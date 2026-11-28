"""只改变第五背景的标尺时，实现五组正交汇所需的独立谱线数。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
SHELL = tuple(2306888100085+j for j in range(5))
HEAVY = tuple(2345666+2*j for j in range(5))
SCALES = (7000,7000,7000,7000,6000)
DIRECTION = (
    (-8173477639,-100000000000,28713094,-76214115,0),
    (-6383341881,-100000000000,474473,-1239660,0),
    (-34593116934,-100000000000,1716671,-4677419,0),
    (100000000000,-100000000000,-536151298,1527688581,0),
    (-90495957423,-100000000000,1717191124,-4953253357,0),
)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def primal():
    solve = load("one-background-specific-line")["solve_layout"]
    layouts = json.loads((ROOT/"ch39-five-split-scale-layouts.json").read_text())
    cases = [("fifth",SHELL[-1:],HEAVY[-1:],[6500],10**12,"12.3454",set()),
             ("five",SHELL,HEAVY,SCALES,10**7,"27.675",{0,323})]
    for site,lower in [(15,"27.59"),(18,"27.46"),(285,"27.74")]:
        cases.append((f"pair-{site}",SHELL,HEAVY,SCALES,10**7,lower,{0,site}))
    for name,shell,heavy,scales,energy,minimum,expected in cases:
        layout = layouts[name]
        base = [arb(k)/11000 for row in layout["base"] for k in row]
        labels = [list(base) for _ in shell]
        independent = set()
        for line,positions in layout["independent"].items():
            i = int(line)
            independent.add(i)
            for row,k in zip(labels,positions,strict=True):
                row[i] = arb(k)/11000
        assert independent == expected
        extra = [(n,i,arb(k)/11000) for n,i,k in layout["extra"]]
        print("独立线、标尺（百万分之一）、共同交汇能标 GeV",sorted(independent),scales,energy,flush=True)
        solve(labels,extra,independent,[minimum]*len(labels),shell,heavy,False,
              reference_gev=energy,scale_millionths=scales)


def inputs():
    gravity = load("two-zero-cutoff-gravity")["gravity_fields"]
    _,lib,backgrounds = gravity(ROOT,HEAVY,tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    paired = load("shared-twelve-node-gauge-gap")["paired_lines"]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    kernel = load("continuous-split-roots")
    sectors = [(paired(ROOT,shell,n),kernel["series"](bg,fields),bg)
               for n,bg in zip(SHELL,backgrounds,strict=True)]
    return sectors,fields,beta,initial,model,kernel


def dual():
    sectors,fields,beta,initial,model,kernel = inputs()
    projection = load("three-background-calibration")["projection"]
    bound,_ = projection(sectors,fields,beta,initial,model,kernel,DIRECTION,None,{0},
                         scale_millionths=SCALES)
    norm = sum(abs(arb(x)/10**11) for row in DIRECTION for x in row[2:4])
    print("仅第 0 线独立时五背景十项规范差的最大绝对值之下界",(bound/norm).str(14),flush=True)
    assert bound/norm > arb("0.61645")


if __name__=="__main__":
    primal()
    dual()
