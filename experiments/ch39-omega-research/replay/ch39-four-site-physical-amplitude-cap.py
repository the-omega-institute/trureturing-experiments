"""任意共同实四站三平面的中心夸克正振幅需要大于一的物理谱权。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def diag(v):
    return acb_mat([[x if i==j else 0 for j in range(len(v))] for i,x in enumerate(v)])


def inner_square(z):
    return abs(z)**2


def column(data):
    return acb_mat([[v] for v in data])


def unit(v):
    r=sum(inner_square(z) for z in v.entries()).sqrt()
    if not r>0:
        return None
    return v*(1/r)


def gap_at_line(y,ckm,mass,cap):
    reps=[ckm*y,y]
    variance=[];means=[];xs=[]
    for line,spectrum in zip(reps,mass,strict=True):
        h=line.entries()
        mod=[inner_square(z) for z in h]
        means.append(sum(p*w for p,w in zip(mod,spectrum,strict=True)))
        v=sum(mod[i]*mod[j]*(spectrum[i]-spectrum[j])**2 for i,j in combinations(range(3),2))
        if not v.is_finite():
            return None
        variance.append(v)
        tangent=[((spectrum[(i+2)%3]-spectrum[(i+1)%3])*h[(i+1)%3]*h[(i+2)%3]).conjugate() for i in range(3)]
        xs.append(column(tangent))
    stable=max(range(2),key=lambda k:float(variance[k].lower()))
    x=unit(xs[stable])
    if x is None:
        return None
    bases=[x,adj(ckm)*x] if stable==0 else [ckm*x,x]
    low=arb(0);high=arb(1)
    for sector,(var,spectrum,s,vx) in enumerate(zip(variance,mass,means,bases,strict=True)):
        p=sum(inner_square(z)*m for z,m in zip(vx.entries(),spectrum,strict=True))
        # s 和 p 均为单位矢量上的期望，直接与真正的谱范围相交。
        lo,hi=spectrum[0],spectrum[2]
        s=s.intersection(lo.union(hi))
        p=p.intersection(lo.union(hi))
        if not p>0 or not s>0:
            return None
        varlo=max(arb(0),var.lower())
        denom=s.upper()
        if sector==0:
            assert cap-s>0
            denom=min(denom,(cap-s).upper())
        l1=varlo/denom**2
        h1=(var+ (spectrum[0]*spectrum[1]*spectrum[2])/p)/s**2
        if not l1.is_finite() or not h1.is_finite():
            return None
        low=max(low,l1.lower())
        high=min(high,h1.upper())
    return low-high


def pencil_cell(base,other,lo,hi,which,ckm,mass,cap,*,line_test=gap_at_line):
    center=(lo+hi)/2
    radius=(hi-lo)/2
    vals,Q=(base+other*center).eig(right=True,algorithm="rump")
    order=sorted(range(3),key=lambda j:vals[j].real.mid())
    vals=[vals[i].real for i in order]
    assert all(left < right for left,right in zip(vals,vals[1:]))
    us=[unit(column([Q[k,j] for k in range(3)])) for j in order]
    U=acb_mat([[v[i,0] for v in us] for i in range(3)])
    R=adj(U)*other*U
    u=arb(0,radius.upper())
    discs=[]
    for i,v in enumerate(vals):
        spread=radius*sum(abs(R[i,j]).upper() for j in range(3) if j!=i)
        discs.append(v+u*R[i,i].real+arb(0,spread.upper()))
    keep=[]
    for i in which:
        rest=[j for j in range(3) if j!=i]
        if any(discs[i].overlaps(discs[j]) for j in rest):
            keep.append(i)
            continue
        # 在中心本征基中把该支的第 i 分量定为一，另两分量解 Schur 方程。
        lam=discs[i]
        H=acb_mat([[(vals[j]-lam if j==k else 0)+u*R[j,k] for k in rest] for j in rest])
        if H.det().contains(0):
            keep.append(i)
            continue
        rhs=column([-u*R[j,i] for j in rest])
        try:
            small=H.solve(rhs)
        except ZeroDivisionError:
            keep.append(i)
            continue
        tangent={j:small[k,0] for k,j in enumerate(rest)}
        v=unit(column([tangent.get(k,acb(1)) for k in range(3)]))
        if v is None:
            keep.append(i)
            continue
        gap=line_test(U*v,ckm,mass,cap)
        if gap is None or not gap>0:
            keep.append(i)
    return keep


def verify(root):
    c=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec=224
    cap=arb("1.22")
    V=acb_mat(c["standard_matrix"](*[arb(v).sin() for v in (".22704",".04275",".00377")],arb("1.139")))
    assert all(abs(z)>0 for z in V.entries())
    mass=[list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    A=adj(V)*diag([z/mass[0][2] for z in mass[0]])*V
    B=diag([z/mass[1][2] for z in mass[1]])
    assert not (A*B-B*A).det().contains(0)
    for base,other in ((A,B),(B,A)):
        stack=[(arb(-1),arb(1),list(range(3)),0)]
        count=0;depth=0
        while stack:
            left,right,which,level=stack.pop()
            keep=pencil_cell(base,other,left,right,which,V,mass,cap)
            count+=1
            depth=max(depth,level)
            if keep:
                assert level<45,(left,right,keep)
                mid=(left+right)/2
                stack.extend(((left,mid,keep,level+1),(mid,right,keep,level+1)))
        print("实射影铅笔一张图：区间数、最大深度、排除的物理谱权上限",count,depth,cap.str(8),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
