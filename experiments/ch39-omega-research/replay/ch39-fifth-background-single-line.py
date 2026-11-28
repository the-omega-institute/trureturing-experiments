"""沿相同编号步长加入第五个背景后，排除任意单条谱线独立分配。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
DIRECTION = (
    -100000000000,-100000000000,372659199,-989192934,
    -100000000000,-100000000000,8230448,-21549532,
    55386128706,-100000000000,-33298862,94544540,
    15973778876,-100000000000,-9680246869,27659628391,
    -78332461623,-100000000000,34615243725,-100000000000,
)


def verify():
    four = run_path(str(ROOT/"ch39-four-background-crossing.py"))
    ctx = four["four_context"]()
    data,fields,_,_,_,kernel,_ = ctx
    load = four["load"]
    _,_,(rows,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(2_345_674,),tight_neutral=True)
    source = load("shared-twelve-node-gauge-gap")
    edges = source["paired_lines"](ROOT,load("native-rh-ray-pvm")["shell_weights"],2_306_888_100_089)
    assert len(edges)==361 and max(z.upper() for _,_,z in rows) < 1
    assert max(z.upper() for _,_,z in rows) < 7*(arb(7)/1000)*min(a.lower() for a,_,_ in edges)
    data.append((edges,kernel["series"](rows,fields),rows))
    groups,constant,norm = load("shared-large-scale-gap")["project_direction"](
        ctx,DIRECTION,kernel_ratio=7)
    lower = four["lower_projection"](ctx,groups,constant,0)
    print("五组精确达到两项引力，只有第 0 线独立分配，十项规范差的最大绝对值下界",
          (lower/norm).str(15))
    assert lower/norm > arb(".5313")
    alone = ([data[-1]],fields,*ctx[2:])
    vector = (-7388046415,-100000000000,1957255449,-5675768631)
    lines,base,size = load("shared-large-scale-gap")["project_direction"](
        alone,vector,kernel_ratio=1)
    assert max(z.upper() for _,_,z in rows) < (arb(7)/1000)*min(a.lower() for a,_,_ in edges)
    # 一个背景下每条线都独立取任意连续位置分布，逐线取最小后相加。
    bound = four["lower_projection"](alone,lines,base,0)/size
    print("单独第五组、全部谱线自由，最大逆耦合差绝对值下界",bound.str(15))
    assert bound > arb(".2359")


if __name__ == "__main__":
    verify()
