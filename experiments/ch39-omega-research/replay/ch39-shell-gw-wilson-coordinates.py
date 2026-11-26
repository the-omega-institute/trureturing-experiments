"""Gauge-invariant coordinates for the matrix potential from two shortest-edge holonomies."""
from itertools import combinations,product
from pathlib import Path
from runpy import run_path

from flint import acb,acb_mat,arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-higgs-low-band.py")))
model,pts,adj,hs = (m[x] for x in ("model","pts","adj","hs"))
count = model["model"]["count"]
identity = m["yuk"]["identity"](3)
height = acb.zeta_zero(43).imag
angles = [(1-arb(d)/26).acos() for d in (3,5,7)]
p,q = (height*(a-angles[0]) for a in angles[1:])


def links():
    short = []
    flat = set()
    for i,j in combinations(range(6),2):
        dot = sum(a*b for a,b in zip(pts[i],pts[j]))
        if dot==25:
            # Orient all three links from the same cyclic orbit.
            start,end = (i,j) if i%2==0 else (j,i)
            short.append(acb_mat(model["link"](pts[start],pts[end],p,q)))
        elif dot!=17:
            flat.add((i,j))
    connected = {0}
    for _ in range(5):
        for i,j in flat:
            if i in connected or j in connected:
                connected.update((i,j))
    assert len(connected)==6 and len(short)==3 and len(flat)==9
    return short


def frame(words):
    X = acb_mat([v.entries() for v in words])
    return adj(X)*X


if __name__=="__main__":
    short = links()
    for pair in combinations(short,2):
        generators = [z for s in pair for z in (s,adj(s))]
        words = [identity]+generators+[a*b for a,b in product(generators,repeat=2)]
        form = frame(words)
        rows = [[form[i,j] for j in range(9)] for i in range(9)]
        assert count(rows,arb(".00743")**2)==0
        assert count(rows,arb(".00745")**2)>0
    all_generators = [z for s in short for z in (s,adj(s))]
    all_words = [identity]+all_generators+[a*b for a,b in product(all_generators,repeat=2)]
    form = frame(all_words)
    rows = [[form[i,j] for j in range(9)] for i in range(9)]
    assert count(rows,arb(".2724")**2)==0
    assert arb(".06098")<q.sin()<arb(".06099")
    print("all two-loop pairs span Mat(3,C) at word length two, gain in (.00743,.00745);",
          "three-loop gain > .2724")
