"""钟差面的共轭半类各带独立复系数时，全部衰减比上的全局竞争。"""
from collections import Counter
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


ETA=arb(1)/4096
CAP=arb("4.48e-9")


def matrix_terms(exact, table, positive):
    mul,conj=exact["product"],exact["conjugate"]
    result=[]
    # 对角元用 |B|² 的比相位，非对角元用 B² 的积相位。
    for conjugate in (True,False):
        counts=Counter()
        for u,v in product(positive,repeat=2):
            w=mul(u,conj(v) if conjugate else v)
            for (k,z),n in table.items():
                counts[k,mul(w,z)]+=n
                counts[tuple(-a for a in k),mul(w,conj(z))]+=n
        result.append(counts)
    return result


def local_kernels(exact, labels, coords, competitor=False):
    mul,conj=exact["product"],exact["conjugate"]
    if competitor:
        u,v=coords
        shifts={(1,0):u,(0,1):v,(1,-1):mul(u,conj(v))}
        shifts.update({tuple(-i for i in key):conj(z) for key,z in tuple(shifts.items())})
    else:
        x,y,d=coords
        shifts={(1,0):x,(0,1):y,(1,-1):d,(-1,0):-x,(0,-1):-y,(-1,1):-d}
    prepared=[]
    for terms in labels:
        row=[]
        for (k,z),n in terms.items():
            if competitor:
                t=mul(z,shifts[k])
                if t==(1,0):
                    row.append((n,arb(0),arb(0)))
                    continue
                angle=exact["unit"](t).arg()
            else:
                angle=exact["unit"](z).arg()+shifts[k]
            s,c=(angle/2).sin_cos()
            row.append((n,s*s,s*c))
        prepared.append(row)
    return prepared


def matrix_at(data,lo,hi):
    fac=(lo/(1+lo)).sqrt().union((hi/(1+hi)).sqrt())
    v=[]
    for j,row in enumerate(data):
        real=imag=arb(0)
        for n,g,p in row:
            if g.is_zero() and p.is_zero():
                real+=n
                continue
            gl=max(arb(0),g.lower())
            gh=g.upper()
            lower=lo/(lo+gh) if lo>0 else arb(0)
            upper=hi/(hi+gl) if gl>0 else arb(1)
            real+=n*lower.union(upper)
            if j:
                if lo+gl>0:
                    bound=(fac*p/(lo.union(hi)+gl.union(gh))).intersection(arb(0,arb(1)/2))
                else:
                    bound=arb(0,arb(1)/2)
                imag+=n*bound
        v.append(acb(real,imag)/2)
    return v


def near_one(exact,counts,clock,dist,bands):
    rival=local_kernels(exact,counts,clock,True)
    hmax=arb(1)/512
    assert (1-arb(".92"))**2/(4*arb(".92"))<hmax
    subtraction=8*(ETA+12*CAP)
    for l,r in bands:
        pending=[(l,r,[(arb(0),hmax)])]
        leaves=nodes=0
        while pending:
            left,right,intervals=pending.pop()
            nodes+=1
            t=left.union(right)
            target=local_kernels(exact,counts,(t,t-dist,dist))
            def check(a,b):
                C=matrix_at(rival,a,b)
                T=matrix_at(target,a,b)
                return (C[0]-T[0]).real>abs(C[1]-T[1])+subtraction
            rest=[]
            while intervals:
                lo,hi=intervals.pop()
                if check(lo,hi):
                    leaves+=1
                elif not check(hi,hi) or not check(lo,lo):
                    rest.append((lo,hi))
                else:
                    assert hi-lo>arb("1e-40"),(lo,hi,left,right)
                    mid=(lo+hi)/2
                    intervals.extend(((lo,mid),(mid,hi)))
            if rest:
                assert right-left>arb("1e-12"),(left,right,rest[:3])
                mid=(left+right)/2
                pending.extend(((left,mid,rest.copy()),(mid,right,rest)))
        print("接近一的衰减比，角段、相位节点与厄米格数",l.str(8),r.str(8),nodes,leaves,flush=True)


def away(exact, table, positive, clock,dist,bands):
    m=600
    keys=((1,0),(0,1),(1,-1))
    bare=[[sum(n*(exact["unit"](z)**q).real for (k,z),n in table.items() if k==key)
           for key in keys] for q in range(1,m+1)]
    B=[sum(exact["unit"](z)**q for z in positive) for q in range(1,m+1)]
    diagonal=[abs(b)**2 for b in B]
    cross=[b*b for b in B]
    assert sum(diagonal[:2])-abs(sum(cross[:2]))>arb(".06")
    x,y=[exact["unit"](z).arg() for z in clock]
    rival=[sum(v*(q*w).cos() for v,w in zip(row,(x,y,x-y),strict=True)) for q,row in enumerate(bare,1)]
    ratio=8*(arb(1)/4)**22/((1-arb(1)/4)*arb(".06"))
    assert (arb(".1")-192*ratio)/(1+ratio)>arb(".09")
    for left,right in bands:
        t=left.union(right)
        raw=[v-sum(w*(q*g).cos() for w,g in zip(row,(t,t-dist,dist),strict=True))-12*CAP
             for q,(v,row) in enumerate(zip(rival,bare),1)]
        assert all(v>arb(".1") for v in raw[:23])
        low=[(v-ETA).lower() for v in raw]
        pending=[(arb(1)/4,arb(".92"))]
        count=0
        while pending:
            l,r=pending.pop()
            pl=ph=arb(1)
            d=arb(0)
            c=acb(0)
            for a,b,g in zip(low,diagonal,cross,strict=True):
                pl*=l
                ph*=r
                w=a*(pl if a>=0 else ph)
                d+=w*b
                c+=w*g
            err=193*8*ph*r/(1-r)
            if d>abs(c)+err:
                count+=1
            else:
                assert r-l>arb("1e-12"), (l,r)
                mid=(l+r)/2
                pending.extend(((l,mid),(mid,r)))
        print("至0.92的衰减比，角段与有限和正定格数",left.str(8),right.str(8),count,flush=True)


def verify(root):
    signed=run_path(str(root/"ch39-signed-coherent-weights.py"))
    prev=run_path(str(root/"ch39-monotone-winding-action.py"))
    ex=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    corners=run_path(str(root/"ch39-curvature-window-density-cap.py"))
    ctx.prec=176
    corners["restrict_corners"]()
    bands=[(arb(l),arb(r)) for l,r in prev["BANDS"]]
    a,b=prev["verify_phases"](root,shell,bands)
    dist=a.union(b)
    table=Counter((k,z) for k,z in ex["enumerate_faces"]().values() if k!=(0,0))
    assert len(table)==32 and set(table.values())=={3}
    assert all(table[k,ex["conjugate"](z)]==n for (k,z),n in table.items())
    p=[z for k,z in table if k==(1,-1) and z[1]>0]
    assert len(p)==2 and all(table[(1,-1),z]==table[(1,-1),ex["conjugate"](z)] for z in p)
    away(ex,table,p,signed["NEW"],dist,bands)
    terms=matrix_terms(ex,table,p)
    near_one(ex,terms,signed["NEW"],dist,bands)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
