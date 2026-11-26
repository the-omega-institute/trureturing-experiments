"""三种有钟面原生振幅的非负相干叠加，全类权单纯形与全衰减区间的竞争证书。"""
from collections import Counter
from functools import cache
from itertools import combinations_with_replacement
from math import isqrt
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, ctx


PAIRS = tuple(combinations_with_replacement(range(3), 2))


def strictly_copositive(entries):
    """在权重单纯形的顶点、边与可能的内部严格极小逐一作整数符号判定。"""
    xx,xy,xd,yy,yd,dd = map(int,entries)
    if min(xx,yy,dd)<=0:
        return False
    for left,middle,right in ((xx,xy,yy),(xx,xd,dd),(yy,yd,dd)):
        if middle<0 and left*right<=middle*middle:
            return False
    # 令第三坐标为 1-x-y。非正定二阶型的最小值在边界也能取得。
    a,b,c = xx-2*xd+dd, xy-xd-yd+dd, yy-2*yd+dd
    p,q = xd-dd, yd-dd
    det=a*c-b*b
    if a<=0 or det<=0:
        return True
    x,y = b*q-c*p, b*p-a*q
    if x<=0 or y<=0 or x+y>=det:
        return True
    return dd*det>c*p*p-2*b*p*q+a*q*q


def effective_terms(exact, source, table, classes):
    mul, conj = exact["product"],exact["conjugate"]
    totals = {}
    for row,(i,j) in enumerate(PAIRS):
        diffs = Counter(mul(u,conj(v)) for u in classes[i] for v in classes[j])
        # 共轭闭合保证 i,j 与 j,i 的整数标签和重数相同。
        assert diffs==Counter({conj(z):w for z,w in diffs.items()})
        for (k,z), n in table.items():
            for d,m in diffs.items():
                label = (k,mul(z,d))
                if label not in totals:
                    totals[label]=[0]*len(PAIRS)
                totals[label][row]+=n*m
    assert len(totals)==10264
    terms=[(k,z,exact["unit"](z).arg()) for k,z in totals]
    weights=np.array(list(totals.values()),dtype=np.int64).T.copy()
    sizes=np.array([len(classes[i])*len(classes[j]) for i,j in PAIRS],dtype=np.int64)
    assert np.array_equal(weights.sum(axis=1),96*sizes)
    shifts=dict(zip(source["KEYS"],(*source["CLOCKS"],mul(source["CLOCKS"][0],conj(source["CLOCKS"][1]))),strict=True))
    rival=[]
    for k,z,_ in terms:
        phase=mul(z,shifts[k])
        rival.append(0 if phase==(1,0) else source["rounded"]((1-exact["unit"](phase).real)/2,True))
    return terms,weights,sizes,np.array(rival,dtype=np.int64)


def cover(source, terms, weights, sizes, rival, distance, bands):
    grid,kernel=source["GRID"],source["KERNEL"]
    hmax=9*grid//16
    assert kernel*hmax+2*grid<2**63 and 96*14**2*(kernel+1)<2**63
    required=np.array([int((int(n)*(arb(1)/256+12*source["DELTA"])*kernel).upper().ceil().unique_fmpz())
                       for n in sizes],dtype=np.int64)
    def evaluation(g,h,upper):
        if h==0:
            return weights[:,g==0].sum(axis=1)*kernel
        denom=h+g
        numer=kernel*h+denom-1 if upper else kernel*h
        return weights@(numer//denom)
    @cache
    def floor_at(h):
        return evaluation(rival,h,False)
    for l,r in bands:
        pending=[(l,r,[(0,hmax)])]
        leaves=nodes=negative=0
        while pending:
            left,right,intervals=pending.pop()
            nodes+=1
            t=left.union(right)
            target=dict(zip(source["KEYS"],(t,t-distance,distance),strict=True))
            vals=np.array([source["rounded"](((ang+target[k])/2).sin()**2,False) for k,_,ang in terms],
                          dtype=np.int64)
            @cache
            def ceil_at(h):
                return evaluation(vals,h,True)
            unresolved=[]
            while intervals:
                lo,hi=intervals.pop()
                lower=floor_at(lo)-ceil_at(hi)-required
                if strictly_copositive(lower):
                    leaves+=1
                    negative+=int(min(lower)<0)
                    continue
                if (hi-lo<=1 or not strictly_copositive(floor_at(hi)-ceil_at(hi)-required)
                    or not strictly_copositive(floor_at(lo)-ceil_at(lo)-required)):
                    unresolved.append((lo,hi))
                else:
                    mid=isqrt(lo*hi) if lo and hi>2*lo else (lo+hi)//2
                    assert lo<mid<hi
                    intervals.extend(((lo,mid),(mid,hi)))
            if unresolved:
                assert right-left>arb("1e-10"), (left,right,unresolved[:8])
                mid=(left+right)/2
                pending.extend(((left,mid,unresolved.copy()),(mid,right,unresolved)))
        print("角段、相位节点、共正末格、含负矩阵条目的末格",l.str(7),r.str(7),nodes,leaves,negative,flush=True)


def verify(root):
    source=run_path(str(root/"ch39-intrinsic-loop-weights.py"))
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
    assert set(table.values())=={3} and len(table)==32
    classes=[[z for k,z in table if k==key] for key in source["KEYS"]]
    assert list(map(len,classes))==[14,14,4]
    assert all(table[k,exact["conjugate"](z)]==n for (k,z),n in table.items())
    tri=lib["terms"](shell)
    assert sum(spoke for k in source["KEYS"] for _,spoke in tri[k])==6
    for k,labels in zip(source["KEYS"],classes,strict=True):
        units={z:exact["unit"](z) for z in labels}
        found=Counter()
        for native,spoke in tri[k]:
            z,=[z for z,u in units.items() if (native-u).contains(0)]
            assert spoke in (0,1)
            found[z]+=1
        assert all(found[z]==table[k,z] for z in labels)
    # 相干首项严格为正。类权非负时合成的首振幅 / 标签总权是三种比值的凸组合。
    first=min(sum(exact["unit"](z).real for z in labels)/len(labels) for labels in classes)
    assert first>arb(".003")
    ratio=(arb(1)/4)**25*4/(3*first**2)
    assert (arb(".2")-192*ratio)/(1+ratio)>arb(".19")
    source["low_rho"](exact,lib,tri,classes[0],distance,bands)
    print("首阶振幅相对标签总权的统一下界",first.lower().str(13),flush=True)
    terms,weights,sizes,rival=effective_terms(exact,source,table,classes)
    cover(source,terms,weights,sizes,rival,distance,bands)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
