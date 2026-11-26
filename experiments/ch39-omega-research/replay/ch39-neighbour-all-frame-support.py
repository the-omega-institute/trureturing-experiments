"""固定中央复四站味空间，在相邻真实谱上覆盖所有共同味基。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def dual_forms(rows, rotation, lib):
    adj = lib["adj"]
    columns=[]
    all_columns=[]
    for row in rows:
        p=adj(rotation)*row*rotation
        all_columns.append(lib["coordinate_t"](p))
        if p[0,0].real+p[1,1].real>0:
            columns.append(lib["coordinate_t"](p))
    assert len(columns)==7
    # 数值余子式只选择整数检验方向，支撑值由原始区间列重新计算。
    vectors=set()
    for triple in combinations(columns,3):
        minors=[arb_mat([[z for j,z in enumerate(row) if j!=i] for row in triple]).det() for i in range(4)]
        normal=[(-1)**i*float(z.mid()) for i,z in enumerate(minors)]
        scale=max(abs(z) for z in normal)
        if scale==0:
            continue
        v=tuple(round(10**6*z/scale) for z in normal)
        vectors.update((v,tuple(-x for x in v)))
    forms=[]
    for direction in sorted(vectors):
        norm=sum(x*x for x in direction)
        v=[arb(x)/arb(norm).sqrt() for x in direction]
        support=sum(max(arb(0),sum(x*y for x,y in zip(v,col,strict=True)).upper()) for col in all_columns)
        forms.append((v,support))
    return forms


def setup(root,offset):
    native=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    geom=run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    ck=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    rows=native["native"](root,2_306_888_100_085+offset)
    ctx.prec=256
    sectors=[dual_forms(sector,geom["four_basis"](),geom) for sector in rows]
    V=acb_mat(ck["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    mass=[list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    up,down=[geom["diag"]([m/row[2] for m in row]) for row in mass]
    A,B=geom["adj"](V)*up*V,down
    assert not (A*B-B*A).det().contains(0)
    return geom,sectors,V,mass,A,B


def lines(base,other,lo,hi,indices,geom):
    adj=geom["adj"]
    center=(lo+hi)/2
    radius=(hi-lo)/2
    vals,Q=(base+other*center).eig(right=True,algorithm="rump")
    order=sorted(range(3),key=lambda j:vals[j].real.mid())
    values=[vals[j].real for j in order]
    assert all(x<y for x,y in zip(values,values[1:]))
    columns=[geom["unit"](acb_mat([[Q[k,j]] for k in range(3)])) for j in order]
    U=acb_mat([[v[i,0] for v in columns] for i in range(3)])
    R=adj(U)*other*U
    delta=arb(0,radius.upper())
    discs=[]
    for i,val in enumerate(values):
        spread=radius*sum(abs(R[i,j]).upper() for j in range(3) if j!=i)
        discs.append(val+delta*R[i,i].real+arb(0,spread.upper()))
    result={}
    for i in indices:
        rest=[j for j in range(3) if j!=i]
        if any(not (discs[j]<discs[i] if j<i else discs[i]<discs[j]) for j in rest):
            continue
        H=acb_mat([[(values[j]-discs[i] if j==k else 0)+delta*R[j,k] for k in rest] for j in rest])
        if H.det().contains(0):
            continue
        try:
            z=H.solve(acb_mat([[-delta*R[j,i]] for j in rest]))
        except ZeroDivisionError:
            continue
        tiny=dict(zip(rest,z.entries(),strict=True))
        v=acb_mat([[tiny.get(j,acb(1))] for j in range(3)])
        length=sum(abs(x)**2 for x in v.entries()).sqrt()
        if length>0:
            result[i]=U*v*(1/length)
    return result


def moments(y,spec):
    prob=[abs(z)**2 for z in y.entries()]
    k=max(range(3),key=lambda j:prob[j].mid())
    prob[k]=prob[k].intersection(1-sum(v for j,v in enumerate(prob) if j!=k))
    s=spec[k]+sum(prob[i]*(spec[i]-spec[k]) for i in range(3) if i!=k)
    weighted=[prob[(i+1)%3]*prob[(i+2)%3]*(spec[(i+1)%3]-spec[(i+2)%3])**2 for i in range(3)]
    return s,prob,sum(weighted)


def phase_cover(constraints):
    stack=[(0,0)]
    count=0
    worst=None
    while stack:
        depth,cell=stack.pop()
        count+=1
        angle=(arb(2*cell+1,1)*arb.pi())/(2**depth)
        sn,cs=angle.sin_cos()
        margin=max((constant+cc*cs+ss*sn).lower() for constant,cc,ss in constraints)
        if margin>arb("1e-12"):
            worst=margin if worst is None else min(worst,margin)
        elif depth<9:
            stack.extend(((depth+1,2*cell),(depth+1,2*cell+1)))
        else:
            return None
    return count,worst


def exclude(y,basis,t,V,mass,forms,geom):
    adj=geom["adj"]
    yrep=[V*y,y]
    stats=[moments(h,spec) for h,spec in zip(yrep,mass,strict=True)]
    remote=1-basis
    s,prob,var=stats[remote]
    if not var>0:
        return None
    norm=var.sqrt()
    w=[None,None]
    w[remote]=norm
    w[basis]=-t*mass[basis][2]/mass[remote][2]*norm
    u,v=arb(".00384"),arb(".0384")
    a2=u*v
    b=(1-a2).sqrt()
    signs=[]
    margins=[]
    for sign in (-1,1):
        chars=[]
        for (mean,_,_),res in zip(stats,w,strict=True):
            chars.extend((mean-sign*res/u,mean+sign*res/v))
        margin=max(z for h in chars for z in (-h.upper(),h.lower()-1))
        if margin>arb("1e-12"):
            margins.append(margin)
        else:
            signs.append(sign)
    if not signs:
        return 0,min(margins)
    h=yrep[remote].entries()
    spec=mass[remote]
    # z 取另一部门在 y 的正交补中的残差，以质量本征基避开与自身的相减。
    shifts=[sum(prob[j]*(spec[i]-spec[j]) for j in range(3) if j!=i) for i in range(3)]
    z=acb_mat([[shift*hj/norm] for shift,hj in zip(shifts,h,strict=True)])
    x=acb_mat([[((spec[(i+2)%3]-spec[(i+1)%3])*h[(i+1)%3]*h[(i+2)%3]).conjugate()/norm] for i in range(3)])
    xrep=[x,adj(V)*x] if remote==0 else [V*x,x]
    zrep=[z,adj(V)*z] if remote==0 else [V*z,z]
    blocks=[]
    for k,(xq,zq,spec) in enumerate(zip(xrep,zrep,mass,strict=True)):
        p=sum(abs(xq[i,0])**2*spec[i] for i in range(3))
        zz=sum(abs(zq[i,0])**2*spec[i] for i in range(3))
        cross=sum(xq[i,0].conjugate()*zq[i,0]*spec[i] for i in range(3))/b
        blocks.append((p,zz,cross))
    count=0
    for sign in signs:
        constraints=[]
        for (p,zz,k),(s,_,_),res,group in zip(blocks,stats,w,forms,strict=True):
            diagonal=(zz-a2*s-(u-v)*sign*res)/(1-a2)
            for (v0,v1,v2,v3),sup in group:
                constant=v0*p+v1*diagonal-sup
                cc=v2*k.real+v3*k.imag
                ss=v3*k.real-v2*k.imag
                if constant.upper()+abs(cc).upper()+abs(ss).upper()>0:
                    constraints.append((constant,cc,ss))
        if not constraints:
            return None
        cert=phase_cover(constraints)
        if cert is None:
            return None
        n,margin=cert
        count+=n
        margins.append(margin)
    return count,min(margins)


def verify(root,offset):
    geom,forms,V,mass,A,B=setup(root,offset)
    for basis,(base,other) in enumerate(((A,B),(B,A))):
        stack=[(arb(-1),arb(1),(0,1,2),0)]
        pairs=phases=depth=0
        smallest=arb(100)
        while stack:
            lo,hi,indices,level=stack.pop()
            pairs+=1
            depth=max(depth,level)
            candidates=lines(base,other,lo,hi,indices,geom)
            keep=[]
            for i in indices:
                proof=exclude(candidates[i],basis,lo.union(hi),V,mass,forms,geom) if i in candidates else None
                if proof is None:
                    keep.append(i)
                else:
                    n,margin=proof
                    phases+=n
                    smallest=min(smallest,margin)
            if keep:
                assert level<32,(offset,basis,lo,hi,keep)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,keep,level+1),(mid,hi,keep,level+1)))
        print("偏移、坐标图、铅笔盒、相位盒、深度、最小分离量",offset,basis,pairs,phases,depth,smallest.str(12),flush=True)


if __name__=="__main__":
    for offset in (-1,1):
        verify(Path(__file__).resolve().parent,offset)
