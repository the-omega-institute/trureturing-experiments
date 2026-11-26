"""All shortest phase-carrying cycles for the negative-energy Dirichlet-to-Neumann pencil."""
from pathlib import Path
import runpy
from flint import arb,acb
m=runpy.run_path(str(Path(__file__).with_name('ch39-arc-relative-determinant.py')))
ell={d:(m['rho']*m['theta'](d)).sqrt() for _,_,d,_ in m['edges']}
classes={(52,504,724),(172,216,724)}
cutoff=sum(ell[d] for d in (152,292,804))+arb('1e-8')
adj=[[] for _ in range(19)]
for i,j,d,link in m['edges']:
    adj[i].append((j,d,link,1))
    adj[j].append((i,d,link.conjugate(),-1))
hits=[]
def walk(path,length,ka,kb,W,ds):
    i=path[-1]
    for j,d,z,sign in adj[i]:
        new=length+ell[d]
        if new>cutoff:continue
        assert new<cutoff
        a,b=ka+sign*(d==724),kb+sign*(d==804)
        if j==path[0]:
            if len(path)>=3 and path[1]<path[-1] and (a or b):
                hits.append((tuple(sorted(ds+[d])),a,b,W*z))
        elif j>path[0] and j not in path:
            walk(path+[j],new,a,b,W*z,ds+[d])
for i in range(19):
    walk([i],arb(0),0,0,acb(1),[])

from collections import defaultdict
buckets=defaultdict(list)
for ds,a,b,W in hits:
    if a<0 or (a==0 and b<0):
        a,b,W=-a,-b,W.conjugate()
    buckets[(ds,a,b)].append(W)
rows=sorted(buckets,key=lambda v:float(sum(ell[d] for d in v[0])))
seen=[]
for ds,a,b in rows:
    length=sum(ell[d] for d in ds)
    spin=sum(buckets[ds,a,b],acb(0))
    assert spin.imag.contains(0)
    factor=-16*ell[ds[0]]*ell[ds[1]]*ell[ds[2]]*spin.real
    print(ds,(a,b),len(buckets[ds,a,b]),length.str(15),factor.str(12))
    if (a,b) not in [r[0] for r in seen]:
        previous=[sum(ell[d] for d in t[0]) for t in rows if t[0]!=ds and
                  float(sum(ell[d] for d in t[0]))<float(length)]
        assert not previous or max(z.upper() for z in previous)<length
        assert not factor.contains(0)
        seen.append(((a,b),length))
assert [r[0] for r in seen]==[(1,0),(1,-1),(0,1)]
for i,row in enumerate(rows[:-1]):
    assert sum(ell[d] for d in row[0])<sum(ell[d] for d in rows[i+1][0])
assert seen[0][1]+2*min(w.lower() for w in ell.values())>cutoff
