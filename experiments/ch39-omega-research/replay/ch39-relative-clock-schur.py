"""Relative clock defect contracts under elimination. Finite 50-site scalar bound."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
clock = run_path(str(root/"ch39-coherent-clock-normalization.py"))
base, B, rounded = [clock[k] for k in ("base","B","rounded")]


def main():
    weights = clock["source"]["tables"]()[0]
    ids = sorted(weights)
    prefix = base["data"]("ch39-tick-separators.json")["prefix"]
    z = [clock["source"]["amb"]["upper"]["height_cell"](prefix[i-1][0]) for i in ids]
    m = sum((weights[i]*v for i,v in zip(ids,z)),B(0))
    target = sum((weights[i]*v*v for i,v in zip(ids,z)),B(0))/(m*m)
    p = [rounded(weights[i]*v/m) for i,v in zip(ids,z)]
    epsilon = B(Q(5,8))
    eps2 = epsilon**2
    cost, trace = B(0), B(1)
    for i,j,terms in clock["clock_pairs"](ids):
        gap = sum(terms,B(0))
        weight = B(2)*clock["root_cell"](p[i]*p[j])
        cost = rounded(cost+weight/gap)
        trace = rounded(trace+weight*eps2/(B(2)*gap+eps2))
    assert Q("5.3999") < cost.lo <= cost.hi < Q("5.4000")
    assert Q("1.6223") < trace.lo <= trace.hi < Q("1.6224")
    bound = eps2*cost*trace/B(2)
    assert bound.hi < Q("1.309")**2 and target.lo > Q("2.410")
    print("For every positive microscopic action with relative four-clock defect at most 5/8,")
    print("the calibrated effective reading is below 2.309 first moments; the alternative exceeds 2.410.")


if __name__=="__main__":
    main()
