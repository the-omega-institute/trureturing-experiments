"""Cyclic shell symmetry and harmonic currents for the lightest Kirchhoff weak doublet."""
from pathlib import Path
import runpy
from flint import acb_mat

weak=runpy.run_path(str(Path(__file__).with_name("ch39-kirchhoff-weak-gauge.py")))
q,pi=weak["q"],weak["pi"]
sites=weak["old"]["sh"]["sites"]
perm=[sites.index((y,z,x)) for x,y,z in sites]
orbits=sorted(set(tuple(sorted({i,perm[i],perm[perm[i]]})) for i in range(19)))
edges={(i,j) for i,j,_ in weak["raw"]}
joined={i for orbit in orbits if len(orbit)==3 and tuple(sorted(orbit[:2])) in edges for i in orbit}
other=[i for i in range(19) if i not in joined and perm[i]!=i]
all_nodes=[[i] for i in range(19)]
def counts(x,large):
    K=acb_mat(19,19)
    for i,j,L in weak["raw"]:
        y=x*(L.upper() if large else L.lower())/weak["speed"]
        assert y<pi
        sn,cs=y.sin_cos()
        K[i,i]+=y*cs/sn
        K[j,j]+=y*cs/sn
        K[i,j]=K[j,i]=-y/sn
    def restricted(partition):
        return weak["inertia"](acb_mat([[sum(K[i,j] for i in a for j in b) for b in partition]
                                         for a in partition]),q(0))
    return [restricted(groups) for groups in (all_nodes,orbits,[[i] for i in other])]
low=counts(q("1.19"),False)
high=counts(q("1.22"),True)
assert low[:2]==[1,1] and high[:2]==[3,1]
assert counts(q("1.19"),True)[0]==1 and counts(q("1.22"),False)[0]==3
assert low[2]==high[2]
print("cyclic vertex orbits",orbits)
print("vertices outside cyclic triangles and fixed point",other)
print("negative counts: complete, invariant, supported on outside",low,high)
# At the fixed hub both modes vanish. Each shortest edge is pointwise
# fixed by a coordinate transposition; the doublet's determinant then flips.
assert sites[14]==(17,17,17)
short_edges=[(i,j) for i,j,d,_ in weak['old']['sh']['edges'] if d==12]
assert len(short_edges)==3
for i,j in short_edges:
    assert any(sites[i][a]==sites[i][b] and sites[j][a]==sites[j][b]
               for a,b in ((0,1),(1,2),(0,2)))
assert all(tuple(sorted((14,i))) in edges for i in range(19) if i!=14)
