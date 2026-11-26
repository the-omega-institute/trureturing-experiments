"""三类实或复相干系数的全域排除与旧竞争点反号见证。"""
from collections import Counter
from functools import cache
from itertools import product
from math import isqrt
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, ctx


NEW = ((89301,-131),(1173,-2))
ETA = arb(1)/4096
DELTA = arb("4.48e-9")


def sine_square(t):
    u=t.sin()
    return u*u


def positive(entries):
    a,b,c,d,e,f=entries
    return a>0 and a*d>b*b and a*d*f+2*b*c*e>a*e*e+d*c*c+f*b*b


def fourier(exact, classes, distance, bands):
    nmax=1800
    pairs=((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))
    A=[[sum((exact["unit"](z)**q).real for z in labels) for labels in classes] for q in range(1,nmax+1)]
    grams=[[row[i]*row[j] for i,j in pairs] for row in A]
    assert positive([sum(row[i]*row[j] for row in A[:3])-arb(".0007")*int(i==j) for i,j in pairs])
    x,y=[exact["unit"](z).arg() for z in NEW]
    reference=[3*sum(a*(q*b).cos() for a,b in zip(row,(x,y,x-y),strict=True)) for q,row in enumerate(A,1)]
    diag=(1,0,0,1,0,1)
    for left,right in bands:
        pending=[(left,right)]
        marginal=arb(".1")
        count=0
        while pending:
            l,r=pending.pop()
            t=l.union(r)
            lower=[v-3*sum(a*(q*z).cos() for a,z in zip(row,(t,t-distance,distance),strict=True))-12*DELTA
                   for q,(v,row) in enumerate(zip(reference,A),1)]
            if all(d>marginal for d in lower[:23]):
                count+=1
                continue
            assert r-l>arb("1e-8")
            mid=(l+r)/2
            pending.extend(((l,mid),(mid,r)))
        # 三项前缀的归一分母至少为 .0007*ρ³||r||²；剩余第24圈起绝对尾界。
        ratio=408*(arb(1)/4)**21/(arb(".0007")*(1-arb(1)/4))
        assert (marginal-192*ratio)/(1+ratio)>arb(".09")
        print("低衰减相位末段",count,flush=True)
        pending=[(left,right,[(arb(1)/4,arb(".98"))])]
        count=nodes=0
        while pending:
            l,r,intervals=pending.pop()
            nodes+=1
            t=l.union(r)
            lows=[(v-3*sum(a*(q*z).cos() for a,z in zip(row,(t,t-distance,distance),strict=True))
                   -12*DELTA-ETA).lower() for q,(v,row) in enumerate(zip(reference,A),1)]
            def cert(p,h):
                pl,ph=arb(1),arb(1)
                sums=[arb(0) for _ in pairs]
                for a,g in zip(lows,grams,strict=True):
                    pl*=p
                    ph*=h
                    # 乘以精确正半定秩一 Gram 矩阵，可向下取标量系数。
                    w=a*(pl if a>=0 else ph)
                    sums=[v+w*z for v,z in zip(sums,g,strict=True)]
                tail=193*408*ph*h/(1-h)
                return positive([x-b*tail for x,b in zip(sums,diag,strict=True)])
            undecided=[]
            while intervals:
                p,h=intervals.pop()
                if cert(p,h):
                    count+=1
                elif not cert(p,p) or not cert(h,h):
                    undecided.append((p,h))
                else:
                    assert h-p>arb("1e-10")
                    mid=(p+h)/2
                    intervals.extend(((p,mid),(mid,h)))
            if undecided:
                assert r-l>arb("1e-8"), (t,undecided[:3])
                mid=(l+r)/2
                pending.extend(((l,mid,undecided.copy()),(mid,r,undecided)))
        print("中段衰减相位节点与正定格",nodes,count,flush=True)


def high(source, coherent, exact, table, classes, distance, bands):
    cfg=dict(source,CLOCKS=NEW)
    terms,weights,sizes,rival=coherent["effective_terms"](exact,cfg,table,classes)
    grid,kernel=source["GRID"],source["KERNEL"]
    hmax=grid//8192
    assert (1-arb(".98"))**2/(4*arb(".98"))<arb(1)/8192
    assert kernel*hmax+2*grid<2**63 and 96*14**2*(kernel+1)<2**63
    required=[int((int(n)*(ETA+12*DELTA)*kernel).upper().ceil().unique_fmpz()) for n in sizes]
    def value(gs,h,upper):
        if h==0:return weights[:,gs==0].sum(axis=1)*kernel
        den=h+gs
        num=kernel*h+den-1 if upper else kernel*h
        return weights@(num//den)
    @cache
    def comp(h,upper):
        return value(rival[upper],h,upper)
    # 对负的 r_i r_j 必须反转该矩阵条目的区间端点。
    shifts=dict(zip(source["KEYS"],(*NEW,exact["product"](NEW[0],exact["conjugate"](NEW[1]))),strict=True))
    comp_low=[]
    for k,z,_ in terms:
        phase=exact["product"](z,shifts[k])
        comp_low.append(0 if phase==(1,0) else source["rounded"]((1-exact["unit"](phase).real)/2,False))
    rival=(rival,np.array(comp_low,dtype=np.int64))
    orthants=[(1,i,j,1,i*j,1) for i,j in product((-1,1),repeat=2)]
    for l,r in bands:
        pending=[(l,r,[(0,hmax)])]
        count=nodes=0
        while pending:
            l,r,intervals=pending.pop()
            nodes+=1
            t=l.union(r)
            angles=dict(zip(source["KEYS"],(t,t-distance,distance),strict=True))
            g=[sine_square((a+angles[k])/2) for k,z,a in terms]
            vals=[np.array([source["rounded"](z,b) for z in g],dtype=np.int64) for b in (True,False)]
            @cache
            def target(h,upper):
                return value(vals[upper],h,upper)
            def cert(lo,hi):
                lower=comp(lo,False)-target(hi,True)
                upper=comp(hi,True)-target(lo,False)
                for signs in orthants:
                    vec=[(int(m) if s>0 else -int(n))-e for m,n,s,e in zip(lower,upper,signs,required,strict=True)]
                    if not coherent["strictly_copositive"](vec):return False
                return True
            unresolved=[]
            while intervals:
                lo,hi=intervals.pop()
                if cert(lo,hi):count+=1
                elif hi-lo<=1 or not cert(hi,hi) or not cert(lo,lo):
                    unresolved.append((lo,hi))
                else:
                    mid=isqrt(lo*hi) if lo and hi>2*lo else (lo+hi)//2
                    assert lo<mid<hi
                    intervals.extend(((lo,mid),(mid,hi)))
            if unresolved:
                assert r-l>arb("1e-10"), (l,r,unresolved[:6])
                mid=(l+r)/2
                pending.extend(((l,mid,unresolved.copy()),(mid,r,unresolved)))
        print("高段衰减相位节点与全正负号格数",nodes,count,flush=True)
    return terms,weights


def flip(root,source,exact,terms,weights,classes,distances):
    h=arb(1)/25000
    rho=(1+2*h)-2*(h*(1+h)).sqrt()
    def kernel(g):return h/(h+g)
    x=arb("-.01046")
    y=x+distances[1]
    theta=(x,y,x-y)
    old=[exact["unit"](z).arg() for z in source["CLOCKS"]]
    refs=(old[0],old[1],old[0]-old[1])
    # r=(1,-1,0)。交叉项以有符号重数准确求和。
    w=weights[0]-2*weights[1]+weights[3]
    gapsum=arb(0)
    for (k,z,a),m in zip(terms,w,strict=True):
        if m:
            i=source["KEYS"].index(k)
            gapsum+=int(m)*(kernel(sine_square((a+refs[i])/2))-kernel(sine_square((a+theta[i])/2)))
    factor=(1+rho)/(2*(1-rho))
    denom=arb(0)
    for i in (0,1):
        for j in (0,1):
            for u,v in product(classes[i],classes[j]):
                phase=exact["product"](u,exact["conjugate"](v))
                denom+=(-1)**(i+j)*(factor*kernel((1-exact["unit"](phase).real)/2)-arb(1)/2)
    assert denom>0
    result=factor*gapsum/denom+arb(0,12*DELTA)
    assert result<0
    physical=run_path(str(root/"ch39-reference-flavor-global-winding.py"))
    physical["physical_all_heights"](root,x,y,DELTA)
    print("旧竞争点反号见证的衰减比、归一作用差",rho.str(15),result.str(12),flush=True)


def verify(root):
    source=run_path(str(root/"ch39-intrinsic-loop-weights.py"))
    coherent=run_path(str(root/"ch39-coherent-intrinsic-weights.py"))
    previous=run_path(str(root/"ch39-monotone-winding-action.py"))
    exact=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    corners=run_path(str(root/"ch39-curvature-window-density-cap.py"))
    ctx.prec=176
    corners["restrict_corners"]()
    bands=[(arb(l),arb(r)) for l,r in previous["BANDS"]]
    distances=previous["verify_phases"](root,shell,bands)
    distance=distances[0].union(distances[1])
    table=Counter((k,z) for k,z in exact["enumerate_faces"]().values() if k!=(0,0))
    assert len(table)==32 and set(table.values())=={3}
    assert all(table[k,exact["conjugate"](z)]==n for (k,z),n in table.items())
    classes=[[z for k,z in table if k==key] for key in source["KEYS"]]
    fourier(exact,classes,distance,bands)
    terms,weights=high(source,coherent,exact,table,classes,distance,bands)
    flip(root,source,exact,terms,weights,classes,distances)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
