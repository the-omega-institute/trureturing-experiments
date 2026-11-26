"""在 1 TeV 质量及 CKM 的联合边际盒上核验物理上型谱峰。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path
from flint import arb, acb, acb_mat, ctx


CENTERS = (".00000615",".00311",".8616",".0000135",".000268",".01401",".22704",".04275",".00377","1.139")
WIDTHS = (".00000014",".00005",".0043",".0000002",".000003",".00009",".00082",".00042",".00008",".023")


def transverse_expectation_floor(spec, own, other, c, d, lam):
    """用前两行的消元关系包围横向正交线的期望，重分量自动给出有利项。"""
    divisor = c * own[0,0] + d * other[0,0] - lam
    if d.contains(0) or other[1,2].contains(0) or divisor.contains(0):
        return spec[0].lower()
    r = other[0,2] / other[1,2]
    diagdiff = c * (own[1,1] - own[0,0]) + d * (other[1,1] - other[0,0])
    denominator = 1 - d * other[1,0] * r / divisor
    if denominator.contains(0):
        return spec[0].lower()
    h_ratio = (r + (diagdiff * r - d * other[0,1]) / divisor) / denominator
    alpha = abs(h_ratio * (spec[2]-spec[0]) / (spec[2]-spec[1]))**2
    a = max(arb(0), alpha.lower())
    return (spec[0] + (spec[1]-spec[0]) * a / (1+a)).lower()


def physical_gap(reps,mass,cap,change,basis,lib,t,lam,sectors):
    low,high=arb(0),arb(1)
    statistics=[]
    for side,(h,spec) in enumerate(zip(reps,mass,strict=True)):
        prob=[abs(x)**2 for x in h.entries()]
        # 利用归一化避免最大分量（常为第三个）的区间半宽进入两次。
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
        own,other=sectors[side]
        c,d=(arb(1),t) if side==basis else (t,arb(1))
        floor=transverse_expectation_floor(spec,own,other,c,d,lam)
        p=p.intersection(floor.union(spec[2]))
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
    var,s,p=statistics[0]
    spec=mass[0]
    self_difference=(var*cap*(2*s-cap)/(cap-s)**2-spec[0]*spec[1]*spec[2]/p)/s**2
    gap=max(gap,self_difference.lower())
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


def checked_cell(lib,ref,ref_other,local,other,change,ref_change,sectors,lo,hi,indices,basis,mass,cap):
    t,h=(lo+hi)/2,(hi-lo)/2
    nominal=ref+t*ref_other
    vals,Q=nominal.eig(right=True,algorithm="rump")
    order=sorted(range(3),key=lambda j:vals[j].real.mid())
    vals=[vals[j].real for j in order]
    assert all(a<b for a,b in zip(vals,vals[1:]))
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
        gap=None if z is None else physical_gap((U*z,H*z) if basis==0 else (H*z,U*z),mass,cap,change,basis,lib,t+u,lam,sectors)
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
    diag,adj=lib["diag"],lib["adj"]
    def mixing(arc):
        return acb_mat(c["standard_matrix"](*[i.sin() for i in arc[:3]],arc[3]))
    @lru_cache(maxsize=512)
    def frame(box,basis):
        params=[arb(mean)+fraction*arb(w)*arb((l+r)/2,(r-l)/2)
                for mean,w,(l,r) in zip(CENTERS,WIDTHS,box,strict=True)]
        mass=[params[:3],params[3:6]]
        assert all(0<row[0] and row[0]<row[1] and row[1]<row[2] for row in mass)
        assert mass[0][2]<cap
        angles=params[6:]
        V=mixing(angles)
        assert all(abs(z)>0 for z in V.entries())
        V0=mixing([z.mid() for z in angles])
        X,Y=[diag([z/row[2] for z in row[:2]]+[arb(1)]) for row in mass]
        x,y=[diag([z[i,i].real.mid() for i in range(3)]) for z in (X,Y)]
        setup=((x,V0*y*adj(V0),X,Y,adj(V),adj(V0)),
               (y,adj(V0)*x*V0,Y,X,V,V0))[basis]
        sectors=((X,V*Y*adj(V)),(Y,adj(V)*X*V))
        return (*setup,sectors),mass
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


def verify_corner(root):
    """在盒的一个角点直接核验共同实四站上的严格正扩张及其上型上界。"""
    ctx.prec = 224
    conventions = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    lib = run_path(str(root / "ch39-four-site-physical-amplitude-cap.py"))
    signs = (1,1,-1,1,1,1,1,-1,-1,-1)
    inputs = [arb(c) + n * arb(w) for c,w,n in zip(CENTERS,WIDTHS,signs,strict=True)]
    up,down = inputs[:3],inputs[3:6]
    V = acb_mat(conventions["standard_matrix"](*[z.sin() for z in inputs[6:9]],inputs[9]))
    A,B = [lib["diag"]([m / spec[2] for m in spec]) for spec in (up,down)]
    vals,Q = (A + arb(".38915") * V * B * lib["adj"](V)).eig(right=True,algorithm="rump")
    order = sorted(range(3),key=lambda j:vals[j].real.mid())
    assert vals[order[2]].real > vals[order[1]].real
    y = lib["unit"](lib["column"]([Q[i,order[2]] for i in range(3)]))
    cap = arb("1.191")
    a_squared = arb(".00090975")
    low,high = arb(0),arb(1)
    for side,(h,spec) in enumerate(((y,up),(lib["adj"](V)*y,down))):
        prob = [abs(z)**2 for z in h.entries()]
        s = sum(pi*m for pi,m in zip(prob,spec,strict=True))
        weighted = [prob[(i+1)%3]*prob[(i+2)%3]*(spec[(i+1)%3]-spec[(i+2)%3])**2 for i in range(3)]
        var = sum(weighted)
        assert var > 0
        p = sum(d*m for d,m in zip(weighted,spec,strict=True))/var
        low = max(low,(var/s**2).upper())
        high = min(high,((var+spec[0]*spec[1]*spec[2]/p)/s**2).lower())
        if side == 0:
            low = max(low,(var/(cap-s)**2).upper())
            det = (cap-spec[0])*(cap-spec[1])*(cap-spec[2])
            high = min(high,((var+det/(cap-p))/(cap-s)**2).lower())
    assert low < a_squared < high,(low,a_squared,high)
    print("四站放宽模型角点的上型上限",cap.str(7),"几何内区间",low.str(14),high.str(14))


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    verify_corner(root)
    verify(root,arb(1),arb("1.18"))
