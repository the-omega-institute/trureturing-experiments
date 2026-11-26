"""在 1 TeV 质量及 CKM 的联合边际盒上核验物理上型谱峰。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path
from flint import arb, acb, acb_mat, ctx


def physical_gap(reps,mass,cap,change,basis,lib,t,lam):
    low,high=arb(0),arb(1)
    statistics=[]
    for side,(h,spec) in enumerate(zip(reps,mass,strict=True)):
        prob=[abs(x)**2 for x in h.entries()]
        # 利用归一化避免成熟分量（常为第三个）的区间半宽进入两次。
        k=max(range(3),key=lambda i:prob[i].mid())
        prob[k]=prob[k].intersection(1-sum(p for j,p in enumerate(prob) if j!=k))
        s=(spec[k]+sum(prob[i]*(spec[i]-spec[k]) for i in range(3) if i!=k))
        s=s.intersection(spec[0].union(spec[-1]))
        weighted=[prob[(i+1)%3]*prob[(i+2)%3]*(spec[(i+1)%3]-spec[(i+2)%3])**2 for i in range(3)]
        var=sum(weighted)
        # 其余两条线的质量及 det/p 可直接用同一对角基中的方差权界住。
        if var>0:
            p=(sum(z*v for z,v in zip(spec,weighted,strict=True))/var).intersection(spec[0].union(spec[2]))
        else:
            p=spec[0].union(spec[2])
        statistics.append((var,s,p))
    remote=1-basis
    ratio=mass[basis][2]/mass[remote][2]
    vl,sl,pl=statistics[basis]
    vr,sr,pr=statistics[remote]
    statistics[basis]=(vl.intersection((ratio*t)**2*vr),sl,pl)
    if not t.contains(0):
        statistics[remote]=(vr.intersection(vl/(ratio*t)**2),sr,pr)
    stable=max(range(2),key=lambda j:statistics[j][0].lower())
    line=reps[stable].entries()
    spec=mass[stable]
    tangent=lib['column']([((spec[(i+2)%3]-spec[(i+1)%3])*line[(i+1)%3]*line[(i+2)%3]).conjugate() for i in range(3)])
    x=lib['unit'](tangent)
    if x is None:
        return None
    other_x=(change if stable==basis else lib['adj'](change))*x
    for side,(var,s,p) in enumerate(statistics):
        spec=mass[side]
        if side!=stable:
            xp=[abs(z)**2 for z in other_x.entries()]
            k=max(range(3),key=lambda i:xp[i].mid())
            cross_p=spec[k]+sum(xp[i]*(spec[i]-spec[k]) for i in range(3) if i!=k)
            p=p.intersection(cross_p)
        statistics[side]=(var,s,p)
        if not p>0 or not s>0:
            return None
        denom=s.upper()
        if side==0:
            assert cap-s>0
            denom=min(denom,(cap-s).upper())
        lower=max(arb(0),var.lower())/denom**2
        upper=(var+spec[0]*spec[1]*spec[2]/p)/s**2
        low=max(low,lower.lower())
        high=min(high,upper.upper())
        if side==0:
            det_complement=(cap-spec[0])*(cap-spec[1])*(cap-spec[2])
            high=min(high,((var+det_complement/(cap-p))/(cap-s)**2).upper())
    gap=low-high
    if t<0:
        # 铅笔本征值精确给出两种归一化期望的差，避免质量矩阵近乎平行时的相减膨胀。
        vl,sl,pl=statistics[basis]
        vr,sr,pr=statistics[remote]
        q=-t*ratio
        for j,vj,sj,pj,si,rho,diff in (
            (remote,vr,sr,pr,sl,q,-mass[basis][2]*lam),
            (basis,vl,sl,pl,sr,1/q,-mass[remote][2]*lam/t)):
            spec=mass[j]
            difference=(vj*diff*(rho*sj+si)/si**2-spec[0]*spec[1]*spec[2]/pj)/sj**2
            gap=max(gap,difference.lower())
    return gap


def checked_cell(lib,ref,ref_other,local,other,change,ref_change,lo,hi,indices,basis,mass,cap):
    t,h=(lo+hi)/2,(hi-lo)/2
    nominal=ref+t*ref_other
    vals,Q=nominal.eig(right=True,algorithm="rump")
    order=sorted(range(3),key=lambda j:vals[j].real.mid())
    vals=[vals[j].real for j in order]
    cols=[lib["unit"](lib["column"]([Q[i,j] for i in range(3)])) for j in order]
    U=acb_mat([[v[i,0] for v in cols] for i in range(3)])
    adj=lib["adj"]
    H0=ref_change*U
    dH=(change-ref_change)*U
    H=H0+dH
    correction0=adj(U)*(local-ref)*U
    cross=adj(H0)*other*dH
    other_mid=lib["diag"]([other[i,i].real.mid() for i in range(3)])
    correction1=adj(H0)*(other-other_mid)*H0+cross+adj(cross)+adj(dH)*other*dH
    u=arb(0,h.upper())
    correction=correction0+t*correction1+u*(adj(H)*other*H)
    discs=[]
    for i,v in enumerate(vals):
        spread=sum(abs(correction[i,j]).upper() for j in range(3) if j!=i)
        discs.append(v+correction[i,i].real+arb(0,spread.upper()))
    keep=[]
    gaps=[]
    merit=1e30
    for i in indices:
        rest=[j for j in range(3) if i!=j]
        separation=min((discs[i].lower()-discs[j].upper() if j<i else discs[j].lower()-discs[i].upper()) for j in rest)
        if not separation>0:
            merit=min(merit,-100+float(separation))
            keep.append(i)
            continue
        lam=discs[i]
        minor=acb_mat([[(vals[j]-lam if j==k else 0)+correction[j,k] for k in rest] for j in rest])
        if minor.det().contains(0):
            merit=min(merit,-50)
            keep.append(i)
            continue
        rhs=lib["column"]([-correction[j,i] for j in rest])
        try:
            small=minor.solve(rhs)
        except ZeroDivisionError:
            merit=min(merit,-50)
            keep.append(i)
            continue
        v=dict(zip(rest,small.entries(),strict=True))
        lam=lam.intersection(vals[i]+correction[i,i].real+sum((correction[i,j]*v[j]).real for j in rest))
        z=lib["unit"](lib["column"]([v.get(k,acb(1)) for k in range(3)]))
        gap=None if z is None else physical_gap((U*z,H*z) if basis==0 else (H*z,U*z),mass,cap,change,basis,lib,t+u,lam)
        merit=min(merit,float(gap.mid()) if gap is not None else -40)
        if gap is not None and gap>arb('1e-8'):
            gaps.append(gap.lower())
        else:
            keep.append(i)
    return keep,gaps,merit


def verify(root,fraction,cap):
    c=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    lib=run_path(str(root/"ch39-four-site-physical-amplitude-cap.py"))
    ctx.prec=224
    centers=(".00000615",".00311",".8616",".0000135",".000268",".01401",".22704",".04275",".00377","1.139")
    widths=(".00000014",".00005",".0043",".0000002",".000003",".00009",".00082",".00042",".00008",".023")
    diag,adj=lib["diag"],lib["adj"]
    def mixing(arc):
        return acb_mat(c["standard_matrix"](*[i.sin() for i in arc[:3]],arc[3]))
    @lru_cache(maxsize=512)
    def frame(box,basis):
        params=[arb(mean)+fraction*arb(w)*arb((l+r)/2,(r-l)/2)
                for mean,w,(l,r) in zip(centers,widths,box,strict=True)]
        mass=[params[:3],params[3:6]]
        angles=params[6:]
        V=mixing(angles)
        assert all(abs(z)>0 for z in V.entries())
        V0=mixing([z.mid() for z in angles])
        X,Y=[diag([z/row[2] for z in row[:2]]+[arb(1)]) for row in mass]
        x,y=[diag([z[i,i].real.mid() for i in range(3)]) for z in (X,Y)]
        setup=((x,V0*y*adj(V0),X,Y,adj(V),adj(V0)),
               (y,adj(V0)*x*V0,Y,X,V,V0))[basis]
        return setup,mass
    results=[]
    for basis in range(2):
        stack=[(((-1.,1.),)*10,-1.,1.,[0,1,2],0)]
        count=0
        depth=0
        minimum=arb(10)
        splits=[0]*11
        def check(box,l,r,which):
            setup,mass=frame(box,basis)
            return checked_cell(lib,*setup,arb(l),arb(r),which,basis,mass,cap)
        while stack:
            box,l,r,which,level=stack.pop()
            keep,gaps,merit=check(box,l,r,which)
            count+=1
            depth=max(depth,level)
            if gaps: minimum=min(minimum,min(gaps))
            if not keep: continue
            assert level<40,(basis,l,r,keep,box)
            mid=(l+r)/2
            axis=10
            if r-l<=1/256:
                best=check(box,mid,mid,keep)[2] if r-l>1/16384 else -1e30
                for i,(lo,hi) in enumerate(box):
                    center=(lo+hi)/2
                    probe=box[:i]+((center,center),)+box[i+1:]
                    score=check(probe,l,r,keep)[2]
                    if score>best:
                        axis=i
                        best=score
            splits[axis]+=1
            if axis==10:
                stack.extend(((box,l,mid,keep,level+1),(box,mid,r,keep,level+1)))
            else:
                lo,hi=box[axis]
                m=(lo+hi)/2
                for interval in ((lo,m),(m,hi)):
                    child=box[:axis]+(interval,)+box[axis+1:]
                    stack.append((child,l,r,keep,level+1))
        results.append((count,depth,minimum,splits))
        print("坐标图",basis,"总区间、深度、最小间隔、各轴二分",results[-1],flush=True)
    print("联合半宽倍数、上型物理谱峰排除上限",fraction.str(6),cap.str(6))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent,arb(1),arb(1))
