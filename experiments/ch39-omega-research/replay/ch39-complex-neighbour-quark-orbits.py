"""固定复物理嵌入，任意无上限非负谱权及任一共同味基下的邻点中心排除。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def rays(root, rotation):
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    result = []
    for index in (2_306_888_100_084,2_306_888_100_086):
        up,_ = native["native"](root,index)
        directions = []
        for p in up:
            trace = (p[0,0]+p[1,1]).real
            if trace > 0:
                t = native["adj"](rotation)*p*rotation
                directions.append((t[0,0].real/trace,t[0,1]/trace))
        assert len(directions) == 7
        result.append(sorted(directions,key=lambda v:v[0].mid()))
    return result


def excluded_by_slice(h, radius_sq, directions):
    # 在指定的第一个对角元上，归一后的凸包切片的极点至多用两根射线。
    if h<directions[0][0] or h>directions[-1][0]:
        return True
    for (h1,z1),(h2,z2) in combinations(directions,2):
        assert h1<h2
        if h<h1 or h>h2:
            continue
        w = ((h-h1)/(h2-h1)).intersection(arb(0).union(1))
        if not abs(z1+(z2-z1)*w)**2 < radius_sq:
            return False
    return True


def statistics(y, spec):
    prob=[abs(z)**2 for z in y.entries()]
    k=max(range(3),key=lambda j:prob[j].mid())
    prob[k]=prob[k].intersection(1-sum(prob[j] for j in range(3) if j!=k))
    s=(spec[k]+sum(prob[j]*(spec[j]-spec[k]) for j in range(3) if j!=k)).intersection(spec[0].union(spec[2]))
    cross=[prob[(j+1)%3]*prob[(j+2)%3]*(spec[(j+1)%3]-spec[(j+2)%3])**2 for j in range(3)]
    variance=sum(cross)
    variance=variance.intersection(arb(0).union((spec[2]-spec[0])**2))
    if variance>0:
        xp=[(v/variance).intersection(arb(0).union(1)) for v in cross]
        j=max(range(3),key=lambda k:xp[k].mid())
        xp[j]=xp[j].intersection(1-sum(x for i,x in enumerate(xp) if i!=j))
        p=(spec[j]+sum(xp[k]*(spec[k]-spec[j]) for k in range(3) if k!=j)).intersection(spec[0].union(spec[2]))
        k_sq=sum(xp[i]*xp[j]*(spec[i]-spec[j])**2 for i,j in combinations(range(3),2))
    else:
        p=spec[0].union(spec[2])
        k_sq=arb(0).union(spec[2]**2)
    magnitude = variance.sqrt() if variance>0 else arb(0).union(variance.upper().sqrt())
    assert magnitude.is_finite()
    return s,magnitude,p,k_sq


def line_excluded(y, spectra, mixing, direction, t, chart):
    u,v = arb(".00384"),arb(".0384")
    data=[statistics(z,spec) for z,spec in zip((mixing*y,y),spectra,strict=True)]
    # 同一铅笔的两个残差模按精确的参数成比例。
    ratio=spectra[chart][2]/spectra[1-chart][2]
    w=data[chart][1].intersection(abs(t)*ratio*data[1-chart][1])
    data[chart]=(data[chart][0],w,data[chart][2],data[chart][3])
    for sign1,sign2 in ((1,1),(1,-1),(-1,1),(-1,-1)):
        if (t>0 and sign1==sign2) or (t<0 and sign1!=sign2):
            continue
        allowed=True
        for sector,((s,magnitude,p,k_sq),sign) in enumerate(zip(data,(sign1,sign2),strict=True)):
            w=sign*magnitude
            plus,minus=s-w/u,s+w/v
            if plus<0 or minus<0:
                allowed=False
                break
            diagonal=(sum(spectra[sector])-s-p-u*v*s-(u-v)*w)/(1-u*v)
            if diagonal<0 or p*diagonal<k_sq/(1-u*v):
                allowed=False
                break
            if sector==0:
                trace=p+diagonal
                if trace>0 and excluded_by_slice(p/trace,k_sq/(1-u*v)/trace**2,direction):
                    allowed=False
                    break
        if allowed:
            return False
    return True


def cell(base, other, lo, hi, remaining, test, lib):
    mid,half=(lo+hi)/2,(hi-lo)/2
    values,vec=(base+other*mid).eig(right=True,algorithm="rump")
    order=sorted(range(3),key=lambda i:values[i].real.mid())
    levels=[values[i].real for i in order]
    assert all(a<b for a,b in zip(levels,levels[1:]))
    cols=[lib["unit"](lib["column"]([vec[j,i] for j in range(3)])) for i in order]
    unitary=acb_mat([[v[i,0] for v in cols] for i in range(3)])
    perturb=lib["adj"](unitary)*other*unitary
    delta=arb(0,half.upper())
    discs=[v+delta*perturb[i,i].real+arb(0,(half*sum(abs(perturb[i,j]) for j in range(3) if j!=i)).upper()) for i,v in enumerate(levels)]
    keep=[]
    for i in remaining:
        rest=[j for j in range(3) if j!=i]
        if any(discs[i].overlaps(discs[j]) for j in rest):
            keep.append(i)
            continue
        minor=acb_mat([[(levels[j]-discs[i] if j==k else 0)+delta*perturb[j,k] for k in rest] for j in rest])
        if minor.det().contains(0):
            keep.append(i)
            continue
        try:
            small=minor.solve(lib["column"]([-delta*perturb[j,i] for j in rest]))
        except ZeroDivisionError:
            keep.append(i)
            continue
        tail=dict(zip(rest,small.entries(),strict=True))
        vec=lib["unit"](lib["column"]([tail.get(j,1) for j in range(3)]))
        if vec is None or not test(unitary*vec,mid+delta):
            keep.append(i)
    return keep


def verify(root):
    reference=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    lib=run_path(str(root/"ch39-four-site-physical-amplitude-cap.py"))
    central=run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    ctx.prec=320
    directions=rays(root,central["four_basis"]())
    ctx.prec=224
    ckm=acb_mat(reference["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    spectra=[list(map(arb,x)) for x in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    a=lib["adj"](ckm)*lib["diag"]([v/spectra[0][2] for v in spectra[0]])*ckm
    b=lib["diag"]([v/spectra[1][2] for v in spectra[1]])
    assert not (a*b-b*a).det().contains(0)
    for offset,direction in zip((-1,1),directions,strict=True):
        for chart,(base,other) in enumerate(((a,b),(b,a))):
            stack=[(arb(-1),arb(1),[0,1,2],0)]
            count,depth=0,0
            while stack:
                lo,hi,lines,level=stack.pop()
                count+=1
                depth=max(depth,level)
                keep=cell(base,other,lo,hi,lines,
                          lambda y,t:line_excluded(y,spectra,ckm,direction,t,chart),lib)
                if keep:
                    assert level<42,(offset,chart,lo,hi,keep)
                    m=(lo+hi)/2
                    stack.extend(((lo,m,keep,level+1),(m,hi,keep,level+1)))
            print("邻点偏移、坐标图、全线覆盖区间数与最深二分",offset,chart,count,depth,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
