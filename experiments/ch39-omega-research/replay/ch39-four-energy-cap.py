"""四组背景同时正交汇时，以对偶投影包住最高共同能标。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
DIRECTION = (
    -357525135052,-5083404107201,198161470,-524349595,0,
    21987573669,-3957353751844,-797322,2048037,0,
    59402369564,1520509847062,-5945682850,776620584,-7392621164,
    2094034197,-16703868289,-1904678,5409610,0,
)


def verify():
    four = run_path(str(ROOT/"ch39-four-background-crossing.py"))
    _,high_witness = four["solve_layout"]("ch39-four-high-energy-layout.json",[".0017"]*4,False,567*10**13)
    assert high_witness == 567*10**13
    ctx = four["four_context"]()
    beta = ctx[2]
    model = four["load"]("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()[0]
    logref = (arb(10**16)/model["mz"]).log()
    components = [10*v for j,v in enumerate(DIRECTION) if j%5!=4]
    colours = [arb(v)/10**10 for v in DIRECTION[4::5]]
    assert all(v <= 0 for v in colours)
    groups,constant,_ = four["load"]("shared-large-scale-gap")["project_direction"](
        ctx,components,kernel_ratio=7,low_colours=colours,common_log=logref)
    lower = four["lower_projection"](ctx,groups,constant,0)
    rate = sum((arb(components[j+2])*(beta[0]-beta[2])+arb(components[j+3])*(beta[1]-beta[2]))/10**11
               for j in range(0,len(components),4))
    assert rate < 0
    cap = 10**16*(lower/rate).exp()
    print("四组全程为正时共同交汇能标的上界（GeV）",cap.str(13))
    assert cap < 579*10**13


if __name__ == "__main__":
    verify()
