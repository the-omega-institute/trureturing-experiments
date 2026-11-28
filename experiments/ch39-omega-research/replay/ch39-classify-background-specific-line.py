"""共同标尺 0.007 时，逐一排除 359 条充当唯一独立分配线的选择。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
EXTRA_PROJECTION = (8835503269,-100000000000,-2267320737,6572777239,
                    -34293662645,-100000000000,2656788433,-7691415892)


def verify():
    high = run_path(str(ROOT/"ch39-shared-large-scale-gap.py"))
    ctx = high["load"]()
    data,_,_,_,_,kernel,shared = ctx
    c = arb(7)/1000
    for edges,_,rows in data:
        assert max(z.upper() for _,_,z in rows) < 7*min((c*a).lower() for a,_,_ in edges)
    def minimum(pair):
        lines = [(c*a,c*span,coef) for a,span,coef in pair]
        return shared["projected_minimum"](kernel["basis"],kernel["S"],lines,arb("1e-7"))
    remaining = list(range(361))
    for y in (shared["DIRECTION"],EXTRA_PROJECTION):
        pairs,constant,norm = high["project_direction"](ctx,y,kernel_ratio=7)
        common = [minimum(pair) for pair in pairs]
        total = constant+sum(common)
        next_candidates, exclusions = [], []
        for i in remaining:
            separate = sum(minimum([line]) for line in pairs[i])
            bound = (total-common[i]+separate)/norm
            if bound > 0:
                exclusions.append(bound.lower())
            else:
                next_candidates.append(i)
        remaining = next_candidates
        print("剩余单线编号，已排除选择的规范差下界",remaining,min(exclusions).str(12),flush=True)
    assert remaining == [0,323]


if __name__ == "__main__":
    verify()
