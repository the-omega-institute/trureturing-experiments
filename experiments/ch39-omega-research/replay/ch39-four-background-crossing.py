"""四个实际背景共用交汇能标时，找出唯一可改变分布的谱线。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
SHELL = tuple(range(2_306_888_100_085,2_306_888_100_089))
HEAVY = tuple(range(2_345_666,2_345_673,2))
DIRECTION = (
    -100000000000,-100000000000,-17106719570,49610242304,
    -100000000000,-100000000000,2444554740,-7048512577,
    -100000000000,-100000000000,34491400403,-100000000000,
    -57280714587,-100000000000,34476817047,-100000000000,
)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def solve_layout(filename,minima,fit_energy,reference_energy):
    previous = load("one-background-specific-line")
    layout = json.loads((ROOT/filename).read_text())
    denominator = layout["denominator"]
    base = [arb(k)/denominator for row in layout["base_numerators"] for k in row]
    labels = [[arb(v)/denominator,*base[1:]] for v in layout["line_zero_bases"]]
    extra = [(n,i,arb(v)/denominator) for n,i,v in layout["extra_numerators"]]
    return previous["solve_layout"](labels,extra,0,minima,SHELL,HEAVY,fit_energy,
                                    reference_gev=reference_energy)


def verify():
    _,energy = solve_layout("ch39-four-background-layout.json",["26.9"]*4,True,10**16)
    assert arb(12391979) < energy < arb(12391981)
    print("四个背景的共同交汇能标 GeV",energy.str(12),flush=True)

    ctx = four_context()
    high = load("shared-large-scale-gap")
    groups,constant,norm = high["project_direction"](ctx,DIRECTION,kernel_ratio=7)
    lower = lower_projection(ctx,groups,constant,323)
    print("第 323 线独立分配时四背景的八项规范差最大绝对值下界",(lower/norm).str(12))
    assert lower/norm > arb("4.119")


def four_context():
    ctx = load("shared-large-scale-gap")["load"]()
    data,fields,_,_,_,kernel,_ = ctx
    _,_,two_more = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,HEAVY[2:],tight_neutral=True)
    source = load("shared-twelve-node-gauge-gap")
    shell = load("native-rh-ray-pvm")["shell_weights"]
    for code,rows in zip(SHELL[2:],two_more,strict=True):
        edges = source["paired_lines"](ROOT,shell,code)
        data.append((edges,kernel["series"](rows,fields),rows))
    c = arb(7)/1000
    for edges,_,rows in data:
        assert max(z.upper() for _,_,z in rows) < 7*c*min(a.lower() for a,_,_ in edges)
        assert max(z.upper() for _,_,z in rows) < 1
    return ctx


def lower_projection(ctx, groups, constant, independent):
    kernel,earlier = ctx[-2:]
    c = arb(7)/1000
    def minimum(lines):
        sites = [(c*a,c*span,coeff) for a,span,coeff in lines]
        return earlier["projected_minimum"](kernel["basis"],kernel["S"],sites,arb("1e-6"))
    lower = constant
    for i,lines in enumerate(groups):
        if i==independent:
            lower += sum(minimum([line]) for line in lines)
        else:
            lower += minimum(lines)
    return lower


if __name__ == "__main__":
    verify()
