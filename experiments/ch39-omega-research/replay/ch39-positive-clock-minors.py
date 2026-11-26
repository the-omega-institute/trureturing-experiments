"""All positive clock stiffnesses in a certified neighbourhood of the 68th zero."""
from itertools import combinations,product
from pathlib import Path
import runpy

import numpy as np
from flint import arb,acb,acb_mat

work=Path(__file__).resolve().parent
edge=runpy.run_path(str(work/"ch39-clock-edge-source.py"))
src=edge["src"]
s=src["shell"]
W,J=src["W"],src["J"]
block=src["block"]
pairs=list(combinations(range(3),2))
entry_scale=2**26
norm_limit=arb("3.4")
ratio=arb(".025")


def integer_box(z):
    a=int((entry_scale*z.lower()).floor().unique_fmpz())
    b=int((entry_scale*z.upper()).ceil().unique_fmpz())
    return a,b


def corners(radius):
    p,q=edge["p"]+arb(0,radius),edge["q"]+arb(0,radius)
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    inv=block(m,J,J).inv()
    mixed=block(m,W,J)
    lift=mixed*inv
    k=block(m,W,W)-lift*block(m,J,W)
    ii,oo=[[W.index(i) for i in ring] for ring in (src["inner"],src["outer"])]
    orig=-block(k,ii,oo)
    vectors={}
    for i,j,d in edge["edges"]:
        v=-block(lift,range(7),(J.index(i),))
        v[W.index(j),0]+=m[i,j].conjugate()
        vectors[i,1+(d==804)]=v
    output=[]
    for choice in product(range(3),repeat=6):
        active=[i for i in range(6) if choice[i]]
        t=acb_mat(orig)
        weight=arb(1)
        if active:
            ids=[J.index(i) for i in active]
            g=block(inv,ids,ids)
            weight=4**len(active)*g.det().real
            assert weight>0
            v=acb_mat([[vectors[i,choice[i]][j,0] for i in active] for j in range(7)])
            t-=block(v,ii,range(len(active)))*g.solve(block(v,oo,range(len(active))).conjugate().transpose())
        herm=-t.conjugate().transpose()*t
        for j in range(3):
            herm[j,j]+=norm_limit**2
        assert all(block(herm,range(n),range(n)).det().real>0 for n in (1,2,3))
        compound=[weight*(t[a,c]*t[b,d]-t[a,d]*t[b,c]) for a,b in pairs for c,d in pairs]
        output.append([integer_box(x) for x in [z.real for z in compound]+[z.imag for z in compound]+[weight]])
    return np.asarray(output,dtype=np.int64).reshape((3,)*6+(19,2))


def test(vertices):
    boxes=vertices.reshape(-1,19,2)
    center=boxes.mean(axis=2)
    data=center[:,:18]/center[:,18:]
    # Choose a real linear functional by an approximate nearest point to the hull.
    x=data.mean(axis=0)
    for _ in range(100):
        diff=data[(data@x).argmin()]-x
        scale=diff@diff
        if scale==0:
            break
        x=x+np.clip(-(x@diff)/scale,0,1)*diff
    w=np.rint(1024*x/max(abs(x))).astype(np.int64)
    # The nuclear norm is at most the sum of the Euclidean row lengths.
    rows=w[:9].reshape(3,3)**2+w[9:].reshape(3,3)**2
    nuclear=sum(arb(int(n)).sqrt() for n in rows.sum(axis=1))
    budget=int((ratio*norm_limit**2*nuclear).upper().ceil().unique_fmpz())
    witness=np.r_[w,-budget]
    assert 19*int(abs(boxes).max())*int(abs(witness).max()) < 2**63
    values=boxes[:,:,0]@np.maximum(witness,0)+boxes[:,:,1]@np.minimum(witness,0)
    return values.min()>0


def bisect(vertices):
    choices=[]
    for j in range(6):
        corner=np.moveaxis(vertices,j,0).reshape(3,-1,19,2).mean(axis=3)
        v=corner[:,:,:18]/corner[:,:,18,None]
        for a,b in pairs:
            size=np.mean(np.linalg.norm(v[a]-v[b],axis=-1))
            choices.append((size,j,a,b))
    _,j,a,b=max(choices)
    local=np.moveaxis(vertices,j,0)
    summed=local[a]+local[b]
    mid=summed//2
    mid[...,1] = -((-summed[...,1])//2)
    child=local.copy()
    child[a]=mid
    other=local.copy()
    other[b]=mid
    return np.moveaxis(child,0,j),np.moveaxis(other,0,j)


if __name__=="__main__":
    pending=[(corners(arb(".000001")),0)]
    count=depth=0
    while pending:
        data,level=pending.pop()
        if test(data):
            count+=1
            depth=max(depth,level)
        else:
            assert level<18
            pending.extend((child,level+1) for child in bisect(data))
    print("all positive stiffnesses: raw ratio > .025; products of triangular cells, max depth",count,depth)
    print("geometric ratio lower bound", (ratio/arb("4.51")).str(12))
