"""Recover both clocks using only simple loops that traverse the same shortest shell edge."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-minimal-triangle-readout.py")))
edges=m["edges"]
theta=m["sh"]["theta"]
def profile(lo,hi):
    # Start along the directed shortest edge. This picks one representative
    # of every undirected simple cycle containing that edge.
    records=[]
    for n in (3,4):
        for tail in permutations([i for i in range(19) if i not in (lo,hi)],n-2):
            path=(lo,hi,*tail)
            links=[tuple(sorted((a,b))) for a,b in zip(path,(*path[1:],lo))]
            if all(e in edges for e in links):
                kind=m["face_type"](path)
                ds=tuple(sorted(edges[e] for e in links))
                records.append((path,kind,ds))
    return records
def check():
    shortest=[(i,j) for (i,j),d in edges.items() if d==12]
    for i,j in shortest:
        R=profile(i,j)
        counts=Counter((len(p),k) for p,k,_ in R)
        assert all(k=="constant" for p,k,_ in R if len(p)==3)
        by=[[v for v in R if v[1]==k] for k in ("s","t","mixed")]
        assert all(by)
        sizes=Counter()
        best={}
        for a in by[0]:
            for b in by[1]:
                for c in by[2]:
                    paths=[a[0],b[0],c[0]]
                    sites=set().union(*paths)
                    size=len(sites)
                    sizes[size]+=1
                    if size==6:
                        support={tuple(sorted(e)) for p in paths for e in zip(p,(*p[1:],p[0]))}
                        ds=tuple(sorted(edges[e] for e in support))
                        best.setdefault(ds,paths)
        assert min(sizes)==6
        if (i,j)==shortest[0]:
            vals={ds:sum(theta(d) for d in ds) for ds in best}
            ds=min(vals,key=lambda d:float(vals[d].mid()))
            assert all(vals[ds]<v for k,v in vals.items() if k!=ds)
            print("loop counts at fixed shortest edge",sorted(counts.items()))
            print("three-loop choices by support size",dict(sorted(sizes.items())))
            print("least total angular length on six vertices",vals[ds].str(13),ds,best[ds])
if __name__=="__main__":
    check()
