"""Compare every degree-three cyclic source to explicit zeta heights at the same stiffness."""
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import runpy
import numpy as np
from flint import acb, acb_poly, arb

edge=runpy.run_path(str(Path(__file__).with_name("ch39-edge-global-hierarchy.py")))
fit, ray, wall=edge["fit"],edge["cyclic"]["ray"],edge["wall"]
adj, I=ray["adjoint"],ray["identity"]
order=13
unit=[acb(0,-2*arb.pi()*k/order).exp() for k in range(order)]
reference=fit["model"](0,0)[0].det().real
assert reference>0
Q=edge["Q"]


def controls(vals):
    poly=acb_poly.interpolate(ray["nodes"],vals)
    return [sum(poly[j]*arb(comb(k,j))/comb(3,j) for j in range(k+1)) for k in range(4)]


def records(p,q):
    m,orbits=fit["model"](p,q)
    det=m.det().real/reference
    result=[]
    for orbit in orbits:
        h, samples=ray["wall_samples"](m,orbit)
        response=adj(h)*m.solve(h)
        score=controls([det*((1-x)*I+4*x*response).det().real for x in ray["nodes"]])
        def bond(i,j):
            return controls([mat[wall.index(i),wall.index(j)] for mat in samples])
        result.append(score+bond(7,6)+bond(11,12))
    return result


def inputs(sample=records):
    grid=[[sample(2*arb.pi()*i/order,2*arb.pi()*j/order) for j in range(order)] for i in range(order)]
    witnesses=[]
    for n in (36,315):
        height=acb.zeta_zero(n).imag
        p,q=(fit["s"][key]*height for key in ("a","b"))
        m=fit["model"](p,q)[0]
        assert all(fit["block"](m,range(k),range(k)).det().real>0 for k in range(1,20))
        witnesses.extend(row[:4] for row in sample(p,q))
    branches=[]
    for orbit in range(4):
        # Along a determinant permutation, each of the six clock-start sites occurs
        # at most once as a row and once as a column, so each frequency lies in [-6,6].
        partial={(i,b):[sum(grid[i][j][orbit][t]*unit[j*b%order] for j in range(order))/order
                        for t in range(12)] for i in range(order) for b in range(-6,7)}
        terms={}
        tail=[arb(0) for _ in range(12)]
        for a,b in product(range(-6,7),repeat=2):
            data=[sum(partial[i,b][t]*unit[i*a%order] for i in range(order))/order for t in range(12)]
            if max(abs(z).upper() for z in data)<arb("0.00002"):
                tail=[v+abs(z) for v,z in zip(tail,data)]
            else:
                terms[a,b]=data
        rounded=edge["rounded"]
        for i,err in enumerate(tail):
            terms[0,0][i]+=acb(arb(0,err.upper()),arb(0,err.upper()))
        branches.append({k:np.array(list(map(rounded,vals)),dtype=np.int64).reshape(3,4,3)
                         for k,vals in terms.items()})
    return branches,np.array([[rounded(z) for z in row] for row in witnesses],dtype=np.int64)


def children(control):
    # Bounds on the four Bernstein controls on either half of a dyadic interval.
    levels=[control]
    while levels[-1].shape[-2]>1:
        a,b=levels[-1][...,:-1,:],levels[-1][...,1:,:]
        z=(a+b)//2
        z[...,2]+=1
        levels.append(z)
    return (np.stack([v[...,0,:] for v in levels],axis=-2),
            np.stack([v[...,-1,:] for v in reversed(levels)],axis=-2))


def stiffness_hull(control):
    low=(control[...,:2]-control[...,2:3]).min(axis=-2)
    high=(control[...,:2]+control[...,2:3]).max(axis=-2)
    mid=(low+high)//2
    return np.concatenate([mid,(high-mid).max(axis=-1,keepdims=True)],axis=-1)


def cover(branch,witnesses,reject=None):
    modes=list(branch)
    @lru_cache(maxsize=2048)
    def at(level,index):
        if level==0:
            return np.stack([branch[k] for k in modes]),witnesses
        c,w=at(level-1,index//2)
        return children(c)[index%2],children(w)[index%2]
    @lru_cache(maxsize=2048)
    def data(level,index):
        c,w=at(level,index)
        wb=stiffness_hull(w)
        return stiffness_hull(c),min(wb[:,0]+wb[:,2])
    @lru_cache(maxsize=32768)
    def ph(level,index,radius):
        return edge["phase"](level,index,radius)
    def check(level,idx):
        coords,inv=np.unique(idx[:,2],return_inverse=True)
        dat=[data(level,int(j)) for j in coords]
        coeff=np.stack([d[0] for d in dat])[inv]
        upper=np.array([d[1] for d in dat])[inv]
        accum=np.zeros((3,len(idx),3),dtype=np.int64)
        for j,(a,b) in enumerate(modes):
            clock=(a*(2*idx[:,0]+1)+b*(2*idx[:,1]+1))%(2<<level)
            ids,undo=np.unique(clock,return_inverse=True)
            phase=np.array([ph(level,int(k),abs(a)+abs(b)) for k in ids])[undo,:,None]
            term=np.stack(edge["times"](coeff[:,j].transpose(2,0,1),phase.transpose(1,0,2)))
            assert np.max(abs(term))<2**51
            accum+=term
        re,im,r=accum
        # At least one explicit witness beats this sector by a uniform relative factor.
        assert np.max(abs(accum))<2**51 and np.max(abs(upper))<2**51
        better=1000*(re[:,0]-r[:,0])>1001*upper
        if reject is not None:
            return better|reject(level,idx,better)
        x=tuple(v[:,1] for v in accum)
        y=tuple(v[:,2] for v in accum)
        xy=edge["times"](x,y)
        return better | (x[1]+x[2]<0) | (xy[1]-xy[2]>0)
    frontier=np.zeros((1,3),dtype=np.int64)
    bits=np.array(list(product((0,1),repeat=3)))
    count=0
    for level in range(15):
        rem=[]
        for idx in np.array_split(frontier,max(1,(len(frontier)+255)//256)):
            done=check(level,idx)
            count+=int(done.sum())
            if not done.all():
                rem.append(idx[~done])
        if not rem:
            return count,level
        frontier=(np.concatenate(rem)[:,None,:]*2+bits).reshape(-1,3)
        assert len(frontier)<3000000,(level,len(frontier))
    raise AssertionError(frontier[0])


if __name__=="__main__":
    branches,witnesses=inputs()
    for i,branch in enumerate(branches):
        print("orbit, retained harmonics, boxes, level",i,len(branch),*cover(branch,witnesses),flush=True)
