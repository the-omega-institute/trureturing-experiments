"""在逐编号独立弱基的中心夸克匹配中检验共同复法向的两个字符份额。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, acb, acb_mat


def strips(box, masses):
    blocks,_,_=box
    inequalities=[]
    for (p,s,_,_,w,det),mass in zip(blocks,masses,strict=True):
        complement=(1-mass[0])*(1-mass[1])*(1-mass[2])
        for x,mean,residual,determinant in ((p,s,w,det),(1-p,1-s,-w,complement)):
            if not x>0 or not mean>0:
                return None
            g=residual/mean
            lower=g**2
            upper=lower+determinant/(x*mean**2)
            inequalities.append((g,lower.lower(),upper.upper()))
    return inequalities


def impossible(rows,lo,hi,amplitude):
    # 当 max(u/v,v/u)<=R 时，|u-v| <= (sqrt(R)-1/sqrt(R))*sqrt(Q)。
    radius=(amplitude*hi.sqrt()).upper()
    lower,upper=-radius,radius
    for slope,left,right in rows:
        required=left.union(right)-lo.union(hi)
        if not slope.contains(0):
            d=required/slope
            lower=max(lower,d.lower())
            upper=min(upper,d.upper())
            if lower>upper:
                return True
        elif not (slope*lower.union(upper)).overlaps(required):
            return True
    return False


def exclude(rows,amplitude,cutoff):
    pending=[(arb(0),arb(1),0)]
    count=0
    while pending:
        lo,hi,level=pending.pop()
        count+=1
        if impossible(rows,lo,hi,amplitude):
            continue
        if level==cutoff:
            return 0
        m=(lo+hi)/2
        pending.extend(((lo,m,level+1),(m,hi,level+1)))
    return count


def verify(root):
    core=run_path(str(root/"ch39-triple-geometry-pencil.py"))
    lib,geom,V,mass,A,B,ratios=core["setup"](root)
    t=arb(".1609")
    Q,D=arb(".00084882"),arb(".013061")
    y=lib["lines"](A,B,t,t,(2,),geom)[2]
    blocks,_,_=core["physical_blocks"](y,0,t,V,mass,geom,lib,ratios)
    for (p,s,_,_,w,det),spec in zip(blocks,mass,strict=True):
        complement=(1-spec[0])*(1-spec[1])*(1-spec[2])
        for x,mean,residual,determinant in ((p,s,w,det),(1-p,1-s,-w,complement)):
            g=residual/mean
            lhs=Q+D*g-g**2
            assert 0<lhs<determinant/(x*mean**2)
    ratio=((D+(D**2+4*Q).sqrt())/2)**2/Q
    assert ratio<arb("1.56")
    print("放宽四站正收缩见证的法向字符份额比",ratio.str(14))
    bounds=(("7","4","1.5595"),("80","18","4"))
    for basis,(a,b) in enumerate(((A,B),(B,A))):
        for branch in range(3):
            bound=arb(bounds[basis][branch])
            amplitude=(bound-1)/bound.sqrt()
            pending=[(arb(-1),arb(1),0)]
            count=depth=geometry=0
            while pending:
                lo,hi,level=pending.pop()
                count+=1
                depth=max(depth,level)
                lines=lib["lines"](a,b,lo,hi,(branch,),geom)
                box=(core["physical_blocks"](lines[branch],basis,lo.union(hi),V,mass,geom,lib,ratios)
                     if branch in lines else None)
                rows=strips(box,mass) if box else None
                cells=exclude(rows,amplitude,min(29,level+17)) if rows else 0
                if cells:
                    geometry+=cells
                    continue
                assert level<25,(basis,branch,float(lo),float(hi),bound)
                mid=(lo+hi)/2
                pending.extend(((lo,mid,level+1),(mid,hi,level+1)))
            print("图、支、份额比下界、铅笔区间、几何区间、最大铅笔深度",
                  basis,branch,bound.str(7),count,geometry,depth,flush=True)


def facet_inequalities(rows, rotation, target, model):
    """七维单位盒的满维四维线性像；每个面法向由三列的精确余子式给出。"""
    adj=model["adj"]
    columns=[]
    for projector in rows:
        block=adj(rotation)*projector*rotation
        if block[0,0].real+block[1,1].real>0:
            columns.append(model["coordinate_t"](block))
    assert len(columns)==7
    assert all(not arb_mat([columns[i] for i in pick]).det().contains(0)
               for pick in combinations(range(7),4))
    k=target[0,1]
    facets=[]
    for pick in combinations(range(7),3):
        n=[(-1)**j*arb_mat([[x for i,x in enumerate(columns[h]) if i!=j]
                             for h in pick]).det() for j in range(4)]
        for sign in (-1,1):
            normal=[sign*z for z in n]
            support=arb(0)
            for h,col in enumerate(columns):
                if h in pick:
                    continue
                value=sum(x*y for x,y in zip(normal,col,strict=True))
                assert value>0 or value<0
                if value>0:
                    support+=value
            n0,n1,n2,n3=normal
            fixed=n0*target[0,0].real+n1*target[1,1].real-support
            facets.append((fixed,n2*k.real+n3*k.imag,n3*k.real-n2*k.imag))
    return facets


def phase_feasible(facets):
    """给定几何与铅笔线，核验单位圆与闭半平面的交集。"""
    if all(a+b<=0 for a,b,c in facets):
        return True
    assert any(a+b>0 for a,b,c in facets)
    for i,(a,b,c) in enumerate(facets):
        norm=b*b+c*c
        assert norm>0
        delta=norm-a*a
        if delta<0:
            continue
        assert delta>0
        for sign in (-1,1):
            shift=sign*delta.sqrt()
            cos=(-a*b-shift*c)/norm
            sin=(-a*c+shift*b)/norm
            others=[d+e*cos+f*sin for j,(d,e,f) in enumerate(facets) if j!=i]
            if all(value<=0 for value in others):
                return True
            assert any(value>0 for value in others)
    return False


def verify_phases(root):
    core=run_path(str(root/"ch39-triple-geometry-pencil.py"))
    lib,geom,V,mass,A,B,ratios=core["setup"](root)
    node=run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    native=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    u,v,_,rotation=node["common_embedding"]()
    weak=native["target_basis"](A,arb(".1")*B)*geom["diag"]([acb(0,1),1,acb(0,-1)])
    targets=[node["extend_target"](geom["adj"](weak)*M*weak*spec[2],u,v)
             for M,spec in zip((A,B),mass,strict=True)]
    assert all(0<block[j,j].real<1 for block in targets for j in (2,3))
    for offset,answer in ((-1,False),(0,True),(1,True)):
        sectors=native["native"](root,2_306_888_100_085+offset)
        facets=[f for rows,target in zip(sectors,targets,strict=True)
                for f in facet_inequalities(rows,rotation,target,geom)]
        found=phase_feasible(facets)
        assert found==answer
        print("固定几何和铅笔线下的实际编号、半平面数、全相位可行性",offset,len(facets),found)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    verify(root)
    verify_phases(root)
