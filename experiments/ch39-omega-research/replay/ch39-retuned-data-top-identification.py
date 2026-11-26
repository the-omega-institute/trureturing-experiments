"""Test all three discrete top orbits against each retuned finite-current data set."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
retune = run_path(str(root/"ch39-right-window-retuned-ground.py"))
fit = run_path(str(root/"ch39-finite-current-width-inverse.py"))
window,study,q = retune["window"],retune["study"],retune["q"]
adj,comp = study["adjoint"],study["comp"]


def projected_pairs(index,windows,pivot):
    paths = [[],[]]
    weak = [window["bounds"](index,upper,windows,pivot,integrals=path)[1]
            for upper,path in zip((False,True),paths,strict=True)]
    initial = [window["flavour"]["before"](index,u) for u in (False,True)]
    start = [(1/(r*a)).union(1/(r*b)) for r,a,b in zip(adj["ratios"],*initial)]
    result = []
    for i in (0,2):
        r,s = adj["ratios"][1],adj["ratios"][i]
        value = start[1]-start[i]
        for (dt,lo,li),(_,hi,ui) in zip(*paths,strict=True):
            state = [x.union(y) for x,y in zip(li,ui)]
            traces = [sum(state[k:k+3]) for k in (3,6,9)]
            # Keep threshold choices correlated within each projected row.
            increments = []
            for flow in (lo,hi):
                corr = sum((b/r-d/s)*g-(c/r-f/s)*t for b,d,g,c,f,t in
                           zip(flow["B"][1],flow["B"][i],state,flow["D"][1],flow["D"][i],traces))
                increments.append(-2*((flow["BETA"][1]/r-flow["BETA"][i]/s)*dt+corr/comp["loop"])/comp["loop"])
            value += increments[0].union(increments[1])
        for stage in windows[-2:]:
            value -= stage[0][1]-stage[0][i]
        result.append(value)
    return weak[0].union(weak[1]),result


def verify():
    initial_time = study["nominal"][2]
    gap = window["graph_gap"]()
    for case,logs in enumerate(retune["log_sets"]):
        _,readings,width = retune["predictions"](logs,initial_time)
        _,*model = fit["invert_observations"](readings,width,index=11431)
        windows = window["window_data"](model,gap)
        exclude_other_orbits(case,windows,q(logs[2]))


def exclude_other_orbits(case,windows,pivot):
    window["verify_speeds"](windows)
    offsets = ((None,"-.00000070","-.0000032"),(".00000070",None,"-.00000265"),
               (".00000315",".0000027",None))
    for index in range(3):
        if index == case:
            continue
        at = pivot+q(offsets[case][index])
        weak,(_,third) = projected_pairs(index,windows,at)
        direction = 1 if index>case else -1
        weak_gap = direction*(study["target_at"](at)-weak)
        pair_gap = -direction*third
        print(case,index,weak_gap.lower().str(10),pair_gap.lower().str(10),flush=True)
        assert weak_gap>arb(2)**-21 and pair_gap>arb(2)**-20
        for upper in (False,True):
            window["bounds"](index,upper,windows,q(33))


if __name__ == "__main__":
    verify()
