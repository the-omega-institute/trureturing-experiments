"""Extend a quarter-turn, degree-three background to all triangular and tetrahedral faces."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
from flint import arb

m=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-flux.py")))
core,base,phases,field=m["main"]()
dot,cross,rot=m["dot"],m["cross"],m["rotate"]
S=sorted(m["shells"][26])
def clique(v):
    return all(26-dot(x,y)<=7 for x,y in combinations(v,2))
edges=[e for e in combinations(S,2) if clique(e)]
faces=[f for f in combinations(S,3) if clique(f)]
tetra=[t for t in combinations(S,4) if clique(t)]
def key(x,y):
    turns=[(x,y),(rot(x),rot(y)),(rot(rot(x)),rot(rot(y)))]
    pos=min(turns)
    neg=min((b,a) for a,b in turns)
    assert pos!=neg
    return (pos,1) if pos<neg else (neg,-1)
orbits=sorted({key(*e)[0] for e in edges})
unknown=[e for e in orbits if e not in field]
assert (len(orbits),len(unknown))==(11,7)
pieces=[]
for v,w,z in faces:
    det=dot(v,cross(w,z))
    den=arb(26).sqrt()*(26+dot(v,w)+dot(w,z)+dot(z,v))
    phi=(arb(det)/den).atan()
    assert abs(phi)<arb.pi()/4
    pieces.append((det,[key(*e) for e in ((v,w),(w,z),(z,v))]))
def turns(det,k):
    k%=4
    if k!=2:
        return k if k<2 else -1
    return 2 if det<0 else (-2 if det>0 else None)
inc=[[(faces.index(t[:i]+t[i+1:]),(-1)**i) for i in range(4)] for t in tetra]
# Bare Berry curvature adds to zero on a tetrahedron since all sites are in the same
# hemisphere and the individual face phases have absolute value strictly less than pi/4.
fixed={e:field[e] for e in orbits if e in field}
accepted=ambiguous=0
for vals in product(range(4),repeat=len(unknown)):
    coeff=fixed|dict(zip(unknown,vals))
    offsets=[turns(det,sum(s*coeff[e] for e,s in row)) for det,row in pieces]
    if None in offsets:
        ambiguous+=1
        continue
    if all(sum(s*offsets[j] for j,s in row)==0 for row in inc):
        accepted+=1
print("fixed core gauge, C3-invariant extensions",4**len(unknown),"with a cut face",ambiguous,
      "admissible with all tetrahedron charges zero",accepted)
assert accepted==3104 and ambiguous==0
