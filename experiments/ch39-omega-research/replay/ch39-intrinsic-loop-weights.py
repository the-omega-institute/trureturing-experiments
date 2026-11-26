"""各钟类及全有钟面原生相位自相关的 Abel 绕数权，随全衰减区间作全局竞争。"""
from collections import Counter
from functools import cache
from math import isqrt
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx
import numpy as np


KEYS = ((1,0),(0,1),(1,-1))
GRID = 2**44
KERNEL = 2**17
DELTA = arb("4.48e-9")
CLOCKS = ((7514,3),(51153,-64))


def spectrum(exact, table, labels):
    mul, conj = exact["product"],exact["conjugate"]
    # 每种相位恰有三张面。供权振幅中共同的因子 3 经归一消去。
    differences = Counter(mul(z,conj(w)) for z in labels for w in labels)
    terms = []
    for key in KEYS:
        counts = Counter()
        for (k,z), n in table.items():
            if k==key:
                for d,w in differences.items():
                    counts[mul(z,d)] += n*w
        assert sum(counts.values()) == len(labels)**2*sum(n for (k,_),n in table.items() if k==key)
        for z,n in counts.items():
            terms.append((key,z,exact["unit"](z).arg(),n))
    return terms


def positive_weights(exact, labels, p, root_mod_p, q):
    """在分裂素数 p 的一个素位上有唯一最小赋值，故每个正整次幂和都非零。"""
    assert p in (7,19) and (root_mod_p**2+867)%p==0 and (2*root_mod_p)%p!=0
    valuations=[]
    for a,b in labels:
        norm=a*a+867*b*b
        e=0
        while norm%p==0:
            e+=1
            norm//=p
        left=(a+b*root_mod_p)%p==0
        right=(a-b*root_mod_p)%p==0
        assert not (left and right) and ((left or right)==(e>0))
        valuations.append(e if left else -e)
    assert valuations.count(min(valuations))==1
    u=[exact["unit"](z) for z in labels]
    assert abs(sum(z**(q+1) for z in u))**2 > 2*abs(sum(z**q for z in u))**2
    print("半衰减比下上升的相邻阶、分裂素数与唯一最小两倍赋值",q,q+1,p,min(valuations),flush=True)


def low_rho(exact, lib, tri, labels, distance, bands):
    """至多四分之一的衰减比用前 25 圈及从 26 圈起的绝对尾界。"""
    x,y = [exact["unit"](z).arg() for z in CLOCKS]
    modes = [lib["coeff"](tri,q,arb(0)) for q in range(1,26)]
    ref = [lib["energy"](q,c,x,y) for q,c in enumerate(modes,1)]
    for left,right in bands:
        pending = [(left,right)]
        while pending:
            l,r=pending.pop()
            t=l.union(r)
            if all(lib["energy"](q,c,t,t-distance)-v > arb(".2")+12*DELTA
                   for q,(c,v) in enumerate(zip(modes,ref),1)):
                continue
            assert r-l>arb("1e-9")
            mid=(l+r)/2
            pending.extend(((l,mid),(mid,r)))
    first=abs(sum(exact["unit"](z) for z in labels))**2
    assert first>0
    tail_ratio=len(labels)**2*(arb(1)/4)**25*4/(3*first)
    assert (arb(".2")-192*tail_ratio)/(1+tail_ratio) > arb(".19")


def rounded(x,upper):
    value = (GRID*x).upper().ceil() if upper else (GRID*x).lower().floor()
    return max(0,min(GRID,int(value.unique_fmpz())))


def full_rho(exact, terms, distance, bands, size):
    # h=(1-rho)^2/(4 rho)。在 rho>=1/4 时，0<h<=9/16。
    mul,conj=exact["product"],exact["conjugate"]
    a,b=CLOCKS
    shifts={k:z for k,z in zip(KEYS,(a,b,mul(a,conj(b))),strict=True)}
    rival=[]
    for k,z,_,_ in terms:
        phase=mul(z,shifts[k])
        rival.append(0 if phase==(1,0) else rounded((1-exact["unit"](phase).real)/2,True))
    rival=np.array(rival,dtype=np.int64)
    weights=np.array([n for _,_,_,n in terms],dtype=np.int64)
    assert int(weights.sum())==96*size**2
    hmax=9*GRID//16
    assert KERNEL*hmax+2*GRID < 2**63 and 96*size**2*(KERNEL+1)<2**63
    # 每个整数下/上界以 2^-17 为单位。H 和 g 的整数单位是 2^-44。
    def count(g,H,upper):
        if H==0:
            return int(weights[g==0].sum())*KERNEL
        denom=H+g
        numer=KERNEL*H+denom-1 if upper else KERNEL*H
        return int(weights@(numer//denom))
    @cache
    def floor_at(H):
        return count(rival,H,False)
    required=(size**2*(arb(1)/256+12*DELTA)*KERNEL).upper().ceil().unique_fmpz()
    required=int(required)
    total_cells=0
    phase_cells=0
    best=96*size**2*KERNEL
    for l,r in bands:
        pending=[(l,r,[(0,hmax)])]
        while pending:
            left,right,intervals=pending.pop()
            phase_cells+=1
            t=left.union(right)
            target={k:shift for k,shift in zip(KEYS,(t,t-distance,distance),strict=True)}
            vals=np.array([rounded(((ang+target[k])/2).sin()**2,False) for k,_,ang,_ in terms],
                          dtype=np.int64)
            @cache
            def ceil_at(H):
                return count(vals,H,True)
            unresolved=[]
            while intervals:
                lo,hi=intervals.pop()
                gap=floor_at(lo)-ceil_at(hi)
                if gap>required:
                    best=min(best,gap)
                    total_cells+=1
                    continue
                if (hi-lo<=1 or floor_at(hi)-ceil_at(hi)<=required
                    or floor_at(lo)-ceil_at(lo)<=required):
                    unresolved.append((lo,hi))
                else:
                    middle=isqrt(lo*hi) if lo and hi>2*lo else (lo+hi)//2
                    assert lo<middle<hi
                    intervals.extend(((lo,middle),(middle,hi)))
            if unresolved:
                assert right-left>arb("1e-9"), (size,left,right,unresolved[:6])
                mid=(left+right)/2
                pending.extend(((left,mid,unresolved.copy()),(mid,right,unresolved)))
    print("供权面组数、合并频率项、相位节点、参数末格、单位核差下界",
          size,len(terms),phase_cells,total_cells,(arb(best)/KERNEL).str(12),flush=True)


def verify(root):
    previous=run_path(str(root/"ch39-monotone-winding-action.py"))
    exact=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    lib=run_path(str(root/"ch39-mincut-global-winding.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    corners=run_path(str(root/"ch39-curvature-window-density-cap.py"))
    ctx.prec=176
    corners["restrict_corners"]()
    bands=[(arb(l),arb(r)) for l,r in previous["BANDS"]]
    d1,d2=previous["verify_phases"](root,shell,bands)
    distance=d1.union(d2)
    table=Counter((k,z) for k,z in exact["enumerate_faces"]().values() if k!=(0,0))
    tri=lib["terms"](shell)
    assert len(table)==32 and set(table.values())=={3}
    assert all(table[k,exact["conjugate"](z)]==n for (k,z),n in table.items())
    for k in KEYS:
        units={z:exact["unit"](z) for key,z in table if key==k}
        found=Counter()
        for native,spoke in tri[k]:
            z,=[v for v,unit in units.items() if (unit-native).contains(0)]
            assert spoke in (0,1)
            found[z]+=1
        assert all(found[z]==table[k,z] for z in found)
    assert sum(spoke for key in KEYS for _,spoke in tri[key])==6
    classes=[[z for k,z in table if k==key] for key in KEYS]
    for labels, witness in zip([*classes,sum(classes,[])],
                               ((7,1,9),(7,6,9),(19,8,1),(7,1,1)),strict=True):
        positive_weights(exact,labels,*witness)
        low_rho(exact,lib,tri,labels,distance,bands)
        terms=spectrum(exact,table,labels)
        full_rho(exact,terms,distance,bands,len(labels))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
