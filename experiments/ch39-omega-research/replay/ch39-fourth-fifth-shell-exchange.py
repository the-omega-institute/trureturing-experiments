"""交换壳谱编号与重谱编号，定位 0.007 标尺上第四、第五背景的可行性变化。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
FOURTH_SHELL = 2_306_888_100_088
FIFTH_SHELL = FOURTH_SHELL+1
HEAVY = tuple(range(2_345_666,2_345_675,2))
SEPARATOR = (-7388046415,-100000000000,1957255449,-5675768631)


def verify():
    four = run_path(str(ROOT/"ch39-four-background-crossing.py"))
    load = four["load"]
    single = load("one-background-specific-line")
    layouts = json.loads((ROOT/"ch39-fourth-fifth-shell-layouts.json").read_text())
    positive = {(FOURTH_SHELL,h) for h in HEAVY} | {(FIFTH_SHELL,HEAVY[2])}
    assert {(v["shell"],v["heavy"]) for v in layouts} == positive
    for layout in layouts:
        labels = [arb(x)/11000 for row in layout["base"] for x in row]
        extra = [(None,i,arb(x)/11000) for i,x in layout["extra"]]
        print("正解的壳编号、重谱编号",layout["shell"],layout["heavy"],flush=True)
        single["solve_layout"]([labels],extra,0,["11.5"],[layout["shell"]],[layout["heavy"]],False,
                               reference_gev=10**12)
    core = load("shared-large-scale-gap")
    ctx = core["load"]()
    fields,kernel = ctx[1],ctx[-2]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT,load("native-rh-ray-pvm")["shell_weights"],FIFTH_SHELL)
    _,_,rows = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,HEAVY,tight_neutral=True)
    floor = arb(7)/1000*min(a.lower() for a,_,_ in edges)
    for number,bg in zip(HEAVY,rows,strict=True):
        if (FIFTH_SHELL,number) in positive:
            continue
        assert max(z.upper() for _,_,z in bg) < min(floor,arb(1))
        single_ctx = ([(edges,kernel["series"](bg,fields),bg)],*ctx[1:])
        groups,base,norm = core["project_direction"](single_ctx,SEPARATOR,kernel_ratio=1)
        lower = four["lower_projection"](single_ctx,groups,base,0)/norm
        print("第五壳编号，重谱编号，最大规范差绝对值下界",number,lower.str(12),flush=True)
        expected = {HEAVY[0]:".8175",HEAVY[1]:".3494",HEAVY[3]:".0885",HEAVY[4]:".2359"}
        assert lower > arb(expected[number])


if __name__ == "__main__":
    verify()
