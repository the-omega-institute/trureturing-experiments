"""Bound every simple open path of the compensating Higgs line on the physical R=26 core."""
from itertools import permutations
from pathlib import Path
from runpy import run_path
from flint import acb,arb


model=run_path(str(Path(__file__).with_name("ch39-shell-gw-symmetric-links.py")))
points=model["points"]
chi=(arb(3).sqrt()+1,acb(1,1))
spin=[]
for i in range(6):
    x,y,z=points[i]
    u,v=arb(26).sqrt()+z,acb(x,y)
    t=chi[0]*u+chi[1].conjugate()*v
    spin.append((u*abs(t)/t,v*abs(t)/t))
edges={}
for i,j in permutations(range(6),2):
    if j!=i^1:
        u,v=spin[i]
        w,z=spin[j]
        edges[i,j]=-2*(u.conjugate()*w+v.conjugate()*z).arg()
angles=[]
for size in range(2,7):
    for path in permutations(range(6),size):
        steps=list(zip(path,path[1:]))
        if all(edge in edges for edge in steps):
            angles.append(abs(sum(edges[edge] for edge in steps)))
lim=max(a.upper() for a in angles)
assert arb(".44787024")<lim<arb(".44787025")
edge=max(abs(a).upper() for a in edges.values())
assert edge<arb(".097057")
def ratio(arc):
    return (arb(".03782")*(arc/2).cos()-arb(".002271")*(arc/2).sin())/arb(".34")
assert ratio(lim)>arb(".1069")
assert ratio(15*edge)>arb(".078")
print("largest simple path phase",lim.str(13),"one doublet ratio bound",ratio(lim).str(12))
print("up to three coupled layers; lower ratio",ratio(15*edge).str(12))
print("extra edge phase necessary to lose the acute-path bound",((arb.pi()/2-lim)/5).str(12))
