"""Eight consecutive RH readings suffice for optimized 41-turn control; seven admit 144."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
control=run_path(str(root/"ch39-two-preparation-pole.py"))
base,B,mul,rounded,root_cell=[control[k] for k in ("base","B","mul","rounded","root_cell")]
prefix=base["data"]("ch39-tick-separators.json")["prefix"]
levels=base["data"]("ch39-minimum-tick-finite-inputs.json")["levels"]
word=base["word"]
cosines=[B(Q(a,word["DEN"]),Q(b,word["DEN"])) for a,b in [word["power"](n)[0] for n in range(42)]]
b41=(B(1)-cosines[41])/B(2)


def node(i):
    if not i:return B(0)
    h=prefix[i-1][0]
    return rounded(B(Q(1<<800,(h+1)**2),Q(1<<800,h*h)))


def quadrature(ids,tol):
    a=[j*j for j in range(1,len(ids))]
    nodes=[node(i) for i in ids]
    weights=[]
    for i,v in enumerate(nodes):
        other=nodes[:i]+nodes[i+1:]
        c=B(Q(1,prod(a)))
        for j,r in enumerate(a,1):
            residue=mul(*(B(Q(-1,r))-x for x in other))/B(prod(1-Q(k,r) for k in a if k!=r))
            t=Q(levels[j],2**288)
            c+=residue*B(t-tol,t+Q(1,2**288)+tol)
        weights.append(rounded(c*mul(*(B(1)+B(k)*v for k in a))/mul(*(v-x for x in other))))
    assert all(w.lo>0 for w in weights)
    return nodes,weights


def boundary(scales,tables,lower):
    p,modes,q2,target=base["physical_inputs"]()
    poles=[B(1)/node(i) for i in (1,2,3)]
    top=poles[2]
    H=B(125**2)*top/q2
    assert H.hi<poles[0].lo and ((B(1)-p)/(B(1)-H*node(1))).hi<1
    expo=mul(*[(top/t)**25 for t in poles[:2]])
    for j,(nodes,wts) in enumerate(tables):
        for lam,mult in zip(modes,base["MODES"]):
            lo=poles[j-1].hi+Q(1,2**36) if j else Q(1,2**30)
            hi=poles[j].lo-Q(1,2**36)
            def fun(y):
                r=sum((rounded(w/(B(1)-B(y)*x)) for x,w in zip(nodes,wts)),B(0))
                return scales*(B(y)-H)-lam+lam*(B(1)-p)*r
            assert fun(lo).hi<0<fun(hi).lo
            for _ in range(40):
                mid=(lo+hi)/2
                value=fun(mid)
                if value.lo>0:hi=mid
                elif value.hi<0:lo=mid
                else:break
            end=hi if lower else lo
            expo=mul(expo,(top/B(end))**mult)
    return rounded(expo**11/target)


def main():
    left_ids=[1,2,3,4,8,9,51,52]
    middle_ids=[0,1,2,3,5,6,15,16]
    A=quadrature(left_ids,Q(1,2**272))
    D=quadrature(middle_ids,Q(1,2**272))
    assert A[1][0].lo>Q(".2166302") and D[1][2].lo>Q(".09") and A[1][2].lo>Q(".06")
    # On the grid away from their roots, the A nodal polynomial is positive
    # and the D nodal polynomial is negative. The ordered pairs prove this
    # for every higher node, including the part not explicitly tabulated.
    assert all(i+1==j for i,j in zip(left_ids[::2],left_ids[1::2]))
    assert middle_ids[:4]==[0,1,2,3] and all(i+1==j for i,j in zip(middle_ids[4::2],middle_ids[5::2]))
    c=rounded(A[1][0]*node(1)*b41)
    residual=boundary(c,(A,D,A),True)
    print("8-row worst-case boundary ratio exceeds",float(residual.lo),flush=True)
    assert residual.lo>1

    ids=[1,2,3,4,11,12,3200]
    table=quadrature(ids,Q(0))
    nodes,weights=table
    actual=B(Q(".000835"))
    residual=boundary(actual,(table,)*3,False)
    print("7-row witness boundary ratio below",float(residual.hi),flush=True)
    assert residual.hi<1
    amplitudes=[root_cell(v*w) for v,w in zip(nodes,weights)]
    a1,a2=amplitudes[:2]
    assert a1.lo>a2.hi and a2.lo>max(v.hi for v in amplitudes[2:])
    loss=B(1)-root_cell((B(1)+cosines[41])/B(2))
    a=sum(amplitudes,B(0))
    assert (a-B(2)*a1-a2*loss).lo>0 and (a-B(2)*a2-a1*loss).lo>0
    sine=[B(0)]+[root_cell((B(1)-v)/B(2)) for v in cosines[1:]]
    for l in range(1,42):
        assert (a1*(sine[l]-sine[l-1])).lo>(a2*sine[1]).hi
    assert (weights[0]*nodes[0]*b41).hi<actual.lo
    assert sum((v*w for v,w in zip(nodes,weights)),B(0)).lo>actual.hi
    parent=run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]
    mean=sum(nodes,B(0))/B(len(nodes))
    for mass,sg in ((69169689551423,1),(69169689551424,-1)):
        assert parent(mass,ids,[dict(zip(ids,weights))],kinetic_scale=actual/mean)==(sg,144)
    print("Positive seven-site countermodel: 69169689551423 < parent mass/GeV < 69169689551424; 144 open.")
    print("No two-preparation allocation of 41 turns gives universal closure on the first seven time rows.")


if __name__=="__main__":
    main()
