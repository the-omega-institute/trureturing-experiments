"""Certify the positive-stiffness mass floor on the entire robust zero-energy clock region."""
from itertools import combinations,product
from functools import lru_cache
from pathlib import Path
import runpy

import numpy as np
from flint import acb,arb,acb_mat

work=Path(__file__).resolve().parent
corn=runpy.run_path(str(work/"ch39-positive-clock-minors.py"))
tri=runpy.run_path(str(work/"ch39-schur-triangulation-density.py"))
edge=corn["edge"]
src=corn["src"]
W,J,block,s=src["W"],src["J"],src["block"],src["shell"]
pairs=corn["pairs"]
pi=arb.pi()
Q=2**21
modes=[(0,0),(1,0),(0,1)]
wide=modes+[(2,0),(1,1),(0,2)]
tests=tri["tests"]+[tri["terms"](tri["loop"](f),tri["y"]) for f in [(6,12,11),(7,12,11)]]


def prepare(z1,z2):
    m=s["matrix"](0,{724:z1,804:z2})
    inv=block(m,J,J).inv()
    lift=block(m,W,J)*inv
    k=block(m,W,W)-lift*block(m,J,W)
    ii,oo=[[W.index(i) for i in ring] for ring in (src["inner"],src["outer"])]
    vertices={}
    for i,j,d in edge["edges"]:
        v=-block(lift,range(7),(J.index(i),))
        v[W.index(j),0]+=m[i,j].conjugate()
        vertices[i,1+(d==804)]=v
    data=[]
    for choice in product(range(3),repeat=6):
        active=[i for i,j in enumerate(choice) if j]
        t=-block(k,ii,oo)
        weight=arb(1)
        if active:
            ids=[J.index(i) for i in active]
            g=block(inv,ids,ids)
            weight=4**len(active)*g.det().real
            v=acb_mat([[vertices[i,choice[i]][j,0] for i in active] for j in range(7)])
            t-=block(v,ii,range(len(active)))*g.solve(block(v,oo,range(len(active))).conjugate().transpose())
        data.append((t,weight))
    return data


def rounded(z):
    a=int((Q*z.lower()).floor().unique_fmpz())
    b=int((Q*z.upper()).ceil().unique_fmpz())
    value=(a+b)//2
    assert max(b-value,value-a)<=2
    return value


def series():
    raw=[]; compound=[]; den=[]
    for (base,weight),(op,_),(oq,_) in zip(prepare(1,1),prepare(-1,1),prepare(1,-1)):
        dx,dy=(base-op)/2,(base-oq)/2
        coeff=[base-dx-dy,dx,dy]
        raw.append([[[rounded(z[i,j].real),rounded(z[i,j].imag)] for i,j in product(range(3),repeat=2)]
                    for z in coeff])
        terms={k:[acb(0)]*9 for k in wide}
        for (i,j),x in zip(modes,coeff):
            for (k,l),y in zip(modes,coeff):
                term=terms[i+k,j+l]
                for v,((a,b),(c,d)) in enumerate(product(pairs,repeat=2)):
                    term[v]+=weight*(x[a,c]*y[b,d]-x[a,d]*y[b,c])
        compound.append([[[rounded(z.real),rounded(z.imag)] for z in terms[k]] for k in wide])
        den.append(rounded(weight))
    return np.array(raw,dtype=np.int64),np.array(compound,dtype=np.int64),np.array(den,dtype=np.int64)


@lru_cache(maxsize=12000)
def unit(level,index):
    z=acb(0,-pi*(index%(2<<level))/(1<<level)).exp()
    return rounded(z.real),rounded(z.imag)


def enclose(data,modes,level,i,j):
    # Integer centers for coefficient and value each miss by at most 2/Q per coordinate.
    angle=int((Q*pi/(1<<level)).upper().ceil().unique_fmpz())
    re=im=error=0
    for k,(a,b) in enumerate(modes):
        c,sine=unit(level,a*(2*i+1)+b*(2*j+1))
        x,y=data[:,k,:,0],data[:,k,:,1]
        re=re+x*c-y*sine
        im=im+x*sine+y*c
        size=abs(x)+abs(y)+4
        error=error+4*Q+2*size+(a+b)*angle*size
    return [(re-error)//Q,-((-(re+error))//Q),(im-error)//Q,-((-(im+error))//Q)]


def cartesian_bounds(data,modes,level,i,j):
    re_lo,re_hi,im_lo,im_hi=enclose(data,modes,level,i,j)
    low=np.concatenate([re_lo,im_lo],axis=1)
    high=np.concatenate([re_hi,im_hi],axis=1)
    return low,high


def proof(boxes,threshold,max_depth):
    todo=[(boxes,0)];count=0
    while todo:
        data,d=todo.pop()
        flat=data.reshape(-1,19,2)
        mid=flat.mean(axis=2);c=mid[:,:18]/mid[:,18:]
        x=c.mean(axis=0)
        for _ in range(36):
            diff=c[(c@x).argmin()]-x
            square=diff@diff
            if square==0:break
            x=x+np.clip(-(x@diff)/square,0,1)*diff
        magnitude=max(abs(x))
        if magnitude==0:return 0
        w=np.rint(1024*x/magnitude).astype(np.int64)
        n=w[:9].reshape(3,3)**2+w[9:].reshape(3,3)**2
        norm=sum(arb(int(v)).sqrt() for v in n.sum(axis=1))
        price=int((norm*threshold).upper().ceil().unique_fmpz())
        witness=np.r_[w,-price]
        assert 19*int(abs(flat).max())*int(abs(witness).max())<2**63
        lower=flat[:,:,0]@np.maximum(witness,0)+flat[:,:,1]@np.minimum(witness,0)
        if lower.min()>0:
            count+=1
        elif d<max_depth:
            todo.extend((child,d+1) for child in corn["bisect"](data))
        else:
            return 0
    return count


def clock_box(t,h,d,l,i,j):
    low,high=cartesian_bounds(t,modes,l,i,j)
    absolute=np.maximum(abs(low),abs(high))
    assert 18*int(absolute.max())**2<2**63
    square=int((absolute**2).sum(axis=1).max())
    norm_squared=arb(square)/Q**2
    lo,hi=cartesian_bounds(h,wide,l,i,j)
    lo=np.column_stack([lo,d-2]);hi=np.column_stack([hi,d+2])
    boxes=np.stack([lo,hi],axis=-1).reshape((3,)*6+(19,2))
    return proof(boxes,norm_squared*arb(".0167"),8 if l<9 else 15)


if __name__=="__main__":
    raw,compound,den=series()
    # Bounds all intermediate integer sums and the squared Frobenius estimate at level >= 6.
    assert max(int(abs(a).max()) for a in (raw,compound,den)) < 12*Q
    assert 100*(12*Q)**2<2**63
    pending=[(0,0,0)];clock_count=source_count=max_level=0
    while pending:
        l,i,j=pending.pop()
        if any(tri["test_value"](test,l,i,j)<0 for test in tests):
            continue
        n=clock_box(raw,compound,den,l,i,j) if l>=6 else 0
        if n:
            clock_count+=1;source_count+=n;max_level=max(max_level,l)
        else:
            assert l<13,(l,i,j)
            pending.extend((l+1,2*i+a,2*j+b) for a,b in product(range(2),repeat=2))
    final=arb(".0167")/arb("4.51")
    assert final>arb(".00316")/arb(".8573")
    print("clock boxes, source cells, max angular depth, normalized ratio floor",
          clock_count,source_count,max_level,final.str(12))
