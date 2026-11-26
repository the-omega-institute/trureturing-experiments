"""Joint optimization over assignments and rotations subject to coplanarity of four arc normals."""
from pathlib import Path
import runpy
from flint import arb,arb_mat

base=runpy.run_path(str(Path(__file__).with_name("ch39-shell-unlabeled-alignment.py")))
S,T,order,geom=[base[k] for k in ("S","target","order","geom")]
den=base["magnitude"]
A=geom["plane_basis"](geom["unit"](S[order.index(0)]))
def constrained(seq):
    sites=[T[j] for j in seq]
    center,first,second=(sites[order.index(j)] for j in (0,7,8))
    assert geom["dots"](center,geom["cross"](first,second))!=0
    B=geom["plane_basis"](geom["unit"](center))
    C=arb_mat([[sum(y[i]*x[j] for x,y in zip(S,sites)) for j in range(3)] for i in range(3)])
    D=B.transpose()*C*A
    scores=[sign*D[2,2]+((D[0,0]+sign*D[1,1])**2+(D[1,0]-sign*D[0,1])**2).sqrt()
            for sign in (1,-1)]
    return max(x.lower() for x in scores).union(max(x.upper() for x in scores))
values={seq:constrained(seq) for seq,_ in base["candidates"]}
top=max(values.values(),key=lambda z:float(z.mid()))
leaders={seq for seq,v in values.items() if not v<top}
assert all(values[seq].overlaps(top) for seq in leaders)
comparison=sorted(((v,seq) for seq,v in values.items() if seq not in leaders),
                  key=lambda row:float(row[0].mid()),reverse=True)
assert top>base["cut"] and top>comparison[0][0]
print("unlabeled coplanar winners",len(leaders),"overlap with geometric winners",len(leaders&set(base["winners"])))
print("minimum excess over global geometric best", (2*(1089-top)/den).str(14))
print("energy above next coplanar candidate", (2*(top-comparison[0][0])/den).str(12))
E0=18-2*1089/den
E2=E0+base["gap"]
margin=(E2.sqrt()-E0.sqrt())/2
print("matching stability threshold: total rms-norm perturbation, per point both shells",
      margin.str(14),(margin/6).str(14))
