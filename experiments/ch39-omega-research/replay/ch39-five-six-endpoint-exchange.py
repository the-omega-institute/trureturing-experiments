"""固定重谱 2345676，交叉使用第五与第六壳的端点和带符号谱权。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
FIFTH, SIXTH = 2306888100089,2306888100090
HEAVY = 2345676


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def primal():
    old = json.loads((ROOT/"ch39-fifth-shell-holdout-layouts.json").read_text())
    native, = [item for item in old if item["heavy"]==HEAVY]
    mixed = json.loads((ROOT/"ch39-five-six-endpoint-layout.json").read_text())
    solve = load("one-background-specific-line")["solve_layout"]
    for layout,weight,c,lower in [(native,FIFTH,7000,"11.6"),(mixed,SIXTH,5500,"11.75")]:
        labels = [arb(k)/11000 for row in layout["base"] for k in row]
        extra = [(None,i,arb(k)/11000) for i,k in layout["extra"]]
        print("正解的端点壳、谱权壳、标尺（百万分之一）",FIFTH,weight,c,flush=True)
        solve([labels],extra,set(),[lower],[FIFTH],[HEAVY],False,reference_gev=10**12,
              scale_millionths=[c],weight_codes=[weight])


def dual():
    full = load("sixth-all-scales")
    _,rows,_,model,initial,beta,fields,_,_ = full["inputs"]()
    kernel = load("continuous-split-roots")
    series = kernel["series"](rows,fields)
    paired = load("shared-twelve-node-gauge-gap")["paired_lines"]
    shell = load("native-rh-ray-pvm")["shell_weights"]
    projection = load("three-background-calibration")["projection"]
    vector1 = (-7364339159,-100000000000,1869761078,-5421825902,0)
    vector2 = (-57642435613,-100000000000,34475176484,-100000000000,0)
    for endpoint,weight,vector,minimum in [(FIFTH,SIXTH,vector1,".16168"),
                                           (SIXTH,FIFTH,vector2,"1.3769"),
                                           (SIXTH,SIXTH,vector2,"1.5220")]:
        edges = paired(ROOT,shell,endpoint,weight_code=weight)
        bound,_ = projection([(edges,series,rows)],fields,beta,initial,model,kernel,[vector],None,set())
        norm = sum(abs(arb(v)/10**11) for v in vector[2:4])
        print("0.007 时端点壳、谱权壳与最大规范差绝对值下界",endpoint,weight,(bound/norm).str(14),flush=True)
        assert bound/norm>arb(minimum)

    swapped = full["inputs"](weight_code=FIFTH)
    assert full["small_scale"](swapped)<arb(".89973")
    parts = json.loads((ROOT/"ch39-sixth-endpoints-fifth-weights-directions.json").read_text())
    assert parts[0][0]==full["START"] and parts[-1][1]==10**6
    assert all(a[1]==b[0] for a,b in zip(parts,parts[1:]))
    results = [full["certify"](part,swapped) for part in parts]
    lower = min(arb(bound).lower() for bound,n in results)
    print("第六壳端点、第五壳谱权的全标尺方向数、闭段数、规范差下界",
          len(parts),sum(n for bound,n in results),lower.str(14))
    assert lower>arb(".03537")


if __name__=="__main__":
    primal()
    dual()
