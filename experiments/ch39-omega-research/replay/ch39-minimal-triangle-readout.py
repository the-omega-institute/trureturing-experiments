"""Minimal neutral triangular Wilson readout and its smallest tetrahedra on the Behrend shell."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-weak-clock-zero-readout.py")))
q=m["q"]
edges=m["edges"]
sh=m["m"]["cur"]["sh"]
clocks={724:"s",804:"t"}
def face_type(face):
    turns=[]
    for i,j in zip(face,(*face[1:],face[0])):
        kind=edges[tuple(sorted((i,j)))]
        if kind in clocks:
            turns.append((clocks[kind],1 if i<j else -1))
    if not turns:
        return "constant"
    if len(turns)==1:
        return turns[0][0]
    assert len(turns)==2 and turns[0][1]!=turns[1][1]
    return "mixed" if turns[0][0]!=turns[1][0] else "constant"

def check():
    tris=[face for face in combinations(range(19),3) if all(e in edges for e in combinations(face,2))]
    kinds={face:face_type(face) for face in tris}
    print("triangle characters",dict(Counter(kinds.values())))
    cells=[]
    for tetra in combinations(range(19),4):
        sides=list(combinations(tetra,2))
        if not all(e in edges for e in sides):
            continue
        types={kinds[f] for f in combinations(tetra,3)}
        if {"s","t","mixed"}<=types:
            lens=tuple(sorted(edges[e] for e in sides))
            cells.append((tetra,lens))
    classes=Counter(v for _,v in cells)
    best=min(classes,key=lambda ds:float(sum(sh["theta"](d) for d in ds).mid()))
    length=sum(sh["theta"](d) for d in best)
    other=min(sum(sh["theta"](d) for d in ds).lower() for ds in classes if ds!=best)
    assert length<other
    selected=[a for a,b in cells if b==best]
    short=[a for a,b in cells if 12 in b]
    print("complete readout tetrahedra",len(cells),"length patterns",len(classes),"with shortest edge",len(short))
    print("minimum total angular length",length.str(12),"gap", (other-length).lower().str(12))
    print("shortest pattern and examples",best,selected[:4],"multiplicity",len(selected))
    if short:
        print("short-edge readout examples",short[:4])
    s,t=m["s"],m["t"]
    axis=m["axis"]
    mixed=m["direct"]-1
    z=2*(1-mixed-(1-s)*(1-t)-axis**2*s*t)/axis
    assert (z**2-16*s*(1-s)*t*(1-t)).contains(0)
    print("719 normalized deficits s, t, mixed",s.str(13),t.str(13),mixed.str(13))
    print("signed invariant sin(2p) sin(2q)",z.str(13))

if __name__=="__main__":
    check()
