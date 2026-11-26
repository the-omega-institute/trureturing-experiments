"""用五个整数时长检验中心恰占五分之一的有序时间表。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
local = run_path(str(root/"ch39-center-free-neighborhood-time.py"))
free = run_path(str(root/"ch39-center-free-three-pulses.py"))


def verify():
    point, scale, eps = local["model"]()
    labels = (11416, 11430, 11431, 11432, 11406)
    counts = (1662932, 8224985, 20000000, 69592639, 519444)
    length = sum(counts)
    assert length == 100000000 and counts[2]*5 == length
    times = [arb(v)/length for v in counts]
    heights, xs, costs = zip(*(point(n) for n in labels))
    gauge = sum(t*x for t,x in zip(times,xs))
    port = gauge + scale*sum(t*c for t,c in zip(times,costs))
    assert abs(gauge) < eps and abs(port) < eps
    print("Normalized gauge and port", (gauge/eps).str(13), (port/eps).str(13))
    ordered, box, study = (free[k] for k in ("ordered", "box", "study"))
    blocks = [study["native"]["plane_response"](h) for h in heights]
    targets = ordered["flavor"]["shell"]["centers"]
    widths = [arb(".000055")]*3+[arb(".006")]
    radius = arb(2)**-28
    outputs = None
    for signs in product((-1,1),repeat=4):
        half = radius/2
        coeff = [c+s*half+arb(0,half) for c,s in zip(box["centers"],signs)]
        steps = {n:ordered["evolution"](block,t,coeff) for n,block,t in zip(labels,blocks,times)}
        readings = ordered["readings"](steps,labels[::-1])
        normalized = [(o-v)/h for o,v,h in zip(readings,targets,widths)]
        assert all(abs(v) < 1 for v in normalized)
        outputs = normalized if outputs is None else [a.union(b) for a,b in zip(outputs,normalized)]
    print("CKM normalized residual bounds", *[(v.lower().str(10),v.upper().str(10)) for v in outputs])


if __name__ == "__main__":
    verify()
