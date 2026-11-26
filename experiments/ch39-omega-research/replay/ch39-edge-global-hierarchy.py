"""All clock angles and positive single-edge strengths, conditional on post-source degree three."""
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import runpy

import numpy as np
from flint import acb, acb_mat, arb

cyclic=runpy.run_path(str(Path(__file__).with_name("ch39-orbit-global-hierarchy.py")))
fit, modes = cyclic["fit"], cyclic["modes"]
block, wall, J = (fit[key] for key in ("block","W","J"))
top_faces=[f for f in fit["src"]["faces"] if f!=(11,12,18)]
top_edges=list(dict.fromkeys((i,j) for f in top_faces for i,j in zip(f,f[1:]+f[:1])))
bottom=(11,12,18)
pairs=list(combinations(range(3),2))
unit=acb_mat([[fit["omega"]**(i*j)/arb(3).sqrt() for j in range(3)] for i in range(3)])
adj=cyclic["ray"]["adjoint"]
Q=2**19


def rounded(z):
    center=np.array([round(float(v.mid())*Q) for v in (z.real,z.imag)],dtype=np.int64)
    radius=int(max((abs(v-arb(int(mid))/Q)*Q).upper().ceil().unique_fmpz()
                   for v,mid in zip((z.real,z.imag),center)))
    return [*center,radius]


def coefficients(orbit):
    values={}
    for i,j in product(range(3),repeat=2):
        m,orbits=fit["model"](2*arb.pi()*i/3,2*arb.pi()*j/3)
        h=orbits[orbit][0]
        hj=block(h,J,(0,))
        c=block(m,J,J).inv()
        proj=block(m,wall,J)*c
        k=block(m,wall,wall)-proj*block(m,J,wall)
        v=block(h,wall,(0,))-proj*hj
        g=(adj(hj)*c*hj)[0,0].real
        assert g>0
        def read(mat):
            t=adj(unit)*block(mat,fit["ii"],fit["oo"])*unit
            return [t[a,b] for a,b in product(range(3),repeat=2)]+[
                mat[wall.index(a),wall.index(b)] for a,b in top_edges]
        correction=-v*adj(v)/g
        base_val=acb(1)
        # All three bottom bonds are independent of both clocks.
        for u,w in zip(bottom,bottom[1:]+bottom[:1]):
            entry=-k[wall.index(u),wall.index(w)]+arb(0).union(arb(1))*correction[wall.index(u),wall.index(w)]
            base_val*=entry
        assert 0<base_val.arg()<arb(".109")
        values[i,j]=[read(-k),read(correction)]
    arr=[]
    for choice in (0,1):
        harmonics=[]
        for a,b in modes:
            factors={(i,j):acb(0,-2*arb.pi()*(i*a+j*b)/3).exp()/9 for i,j in values}
            row=[rounded(sum(z*values[key][choice][t] for key,z in factors.items())) for t in range(9+len(top_edges))]
            harmonics.append(row)
        arr.append(harmonics)
    return np.array(arr,dtype=np.int64)


# Box arithmetic: (a,b,r) represents a/Q+i b/Q with error <=r/Q on each component.
# Integer operations are exact. A bound is checked before each fixed-point multiplication.
def times(z,w):
    a,b,r=z
    c,d,s=w
    A=max(int(np.max(np.abs(t))) for t in z)
    B=max(int(np.max(np.abs(t))) for t in w)
    assert 8*A*B+2*Q<2**63
    error=(np.abs(a)+np.abs(b))*s+(np.abs(c)+np.abs(d))*r+2*r*s
    return (a*c-b*d)//Q,(a*d+b*c)//Q,(error+2*Q-1)//Q


@lru_cache(maxsize=4096)
def phase(level,index,radius):
    turn=(arb(index)+arb(0,radius))/(2<<level)
    return rounded(acb(0,2*arb.pi()*turn).exp())


def field(data,lev,idx):
    # Equal dyadic subdivisions of the two circle coordinates and the compact strength.
    n=1<<lev
    centers=2*idx+1
    result=[]
    for h,(a,b) in enumerate(modes):
        clock=(a*centers[:,0]+b*centers[:,1])%(2*n)
        unique,reverse=np.unique(clock,return_inverse=True)
        ph=np.array([phase(lev,int(j),abs(a)+abs(b)) for j in unique])[reverse,:,None]
        assert n<=Q//2
        stiffness=[centers[:,2,None]*(Q//(2*n)),0,Q//(2*n)]
        first=data[0,h].T
        second=times(data[1,h].T,stiffness)
        coeff=[z+w for z,w in zip(first,second)]
        result.append(times(coeff,[ph[:,j] for j in range(3)]))
    # Coefficients are small: also guard the unscaled integer sum.
    assert all(np.max(np.abs(comp))<2**56 for term in result for comp in term)
    return [sum(term[j] for term in result) for j in range(3)]


def modulus_squared(z,upper):
    # Keep squared lengths independently scaled by Q^2.
    a,b,r=z
    x,y=(np.abs(a)+r,np.abs(b)+r) if upper else (
        np.maximum(0,np.abs(a)-r),np.maximum(0,np.abs(b)-r))
    assert max(np.max(x),np.max(y))<2**26
    return x*x+y*y


def verify_boxes(data,lev,idx):
    a,b,r=field(data,lev,idx)
    passed=np.zeros(len(idx),dtype=bool)
    for f in top_faces:
        # Degree three forces each of six non-bottom faces into (pi-.109,pi).
        loop=(Q,0,0)
        for i,j in zip(f,f[1:]+f[:1]):
            t=9+top_edges.index((i,j))
            loop=times(loop,(a[:,t],b[:,t],r[:,t]))
        re,im,err=loop
        assert max(int(np.max(np.abs(v))) for v in loop)<2**57
        passed|=(im+err<0)|(re-err>0)|(9*(im-err)+(re-err)>0)
    needed=~passed
    if not np.any(needed):
        return passed
    t=[v[needed,:9].reshape(-1,3,3) for v in (a,b,r)]
    frob=modulus_squared(t,True).sum(axis=(1,2))
    assert max(frob)<2**54
    # For the exterior square, any column or row norm is <= s1*s2.
    exterior=[]
    for (i,j),(k,l) in product(pairs,repeat=2):
        def entry(h,v):
            return [x[:,h,v] for x in t]
        z=times(entry(i,k),entry(j,l))
        w=times(entry(i,l),entry(j,k))
        exterior.append(modulus_squared((z[0]-w[0],z[1]-w[1],z[2]+w[2]),False))
    squares=np.stack(exterior,axis=1).reshape(-1,3,3)
    norm2=np.maximum(squares.sum(axis=1).max(axis=1),squares.sum(axis=2).max(axis=1))
    # ceil(0.074*frob/Q) in fixed-point units; all operands are nonnegative.
    assert max(frob)<2**54 and int(max(frob))*37<2**62
    budget=(37*frob+500*Q-1)//(500*Q)
    assert max(budget)<2**30
    passed[needed]=norm2>budget*budget
    return passed


def cover(orbit):
    data=coefficients(orbit)
    frontier=np.zeros((1,3),dtype=np.int64)
    children=np.array(list(product((0,1),repeat=3)))
    count=0
    for level in range(17):
        remaining=[]
        for idx in np.array_split(frontier,max(1,(len(frontier)+511)//512)):
            done=verify_boxes(data,level,idx)
            count+=int(done.sum())
            if not done.all():
                remaining.append(idx[~done])
        if not remaining:
            return count,level
        frontier=(np.concatenate(remaining)[:,None,:]*2+children).reshape(-1,3)
        assert len(frontier)<10000000,(level,orbit)
    raise AssertionError((orbit,frontier[0]))


if __name__=="__main__":
    metric=runpy.run_path(str(Path(__file__).with_name("ch39-broken-edge-maxwell-metric.py")))
    kappa=(arb(793)/473)*metric["oscillation"].upper().exp()
    distortion=kappa*(1+2*kappa)/(kappa+2)
    assert distortion<arb("3.707") and arb(".074")/distortion>arb(".01996")
    assert arb(1)/9>arb(".109").tan()
    print("metric condition product and normalized ratio",distortion.upper().str(12),
          (arb(".074")/distortion).lower().str(12))
    for orbit in range(4):
        print("edge orbit, number of boxes, maximum level",orbit,*cover(orbit),flush=True)
