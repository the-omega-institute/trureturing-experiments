"""第21、22枚实际零点的全部非负双边权不会给出行列式三角主值陈数。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def top_frame(matrix):
    vals, vectors=matrix.eig(right=True,algorithm="rump")
    ids=sorted(range(len(vals)),key=lambda j:-float(vals[j].real))
    assert all(z.imag.contains(0) for z in vals)
    assert all(vals[i].real>vals[j].real for i,j in zip(ids[:3],ids[1:4]))
    cols=[]
    for j in ids[:3]:
        v=[vectors[i,j] for i in range(19)]
        norm=sum(abs(z)**2 for z in v).sqrt()
        cols.append([z/norm for z in v])
    gap=min((vals[i].real-vals[j].real).lower() for i in ids[:3] for j in ids[3:])
    return acb_mat(cols).transpose(),gap


def ingredients(shell,n):
    H=shell["laplacian"](-3,acb.zeta_zero(n).imag)
    B=acb_mat(19,2)
    for k,(i,j) in enumerate(((0,2),(5,14))):
        B[i,k]=1
        B[j,k]=H[j,i]
        assert abs(H[j,i]).contains(1)
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    K=(H-B*adj(B)+eye).inv()
    W=K*B
    return K,W,adj(B)*W


def s_matrix(G,u,v):
    t=acb_mat([[u,0],[0,v]])
    d=acb_mat([[(1-u)/8,0],[0,(1-v)/8]])
    return t*(d+G*t).inv()


def certify(shell,n,anchor):
    K,W,G=ingredients(shell,n)
    ww=adj(W)*W
    @cache
    def S(i,j):
        return s_matrix(G,arb(i)/65536,arb(j)/65536)
    stack=[(0,0,65536)]
    visited=0
    worst=arb(0)
    while stack:
        i,j,d=stack.pop()
        visited+=1
        small,large=S(i,j),S(i+d,j+d)
        shift=((large-small)*ww).trace().real.upper()/2
        center=K-W*((small+large)/2)*adj(W)
        frame,gap=top_frame(center)
        keep=False
        if gap>2*shift:
            err=shift/(gap-shift)
            f=abs((adj(anchor)*frame).det())
            ds=(1-f.lower()**2).max(arb(0)).sqrt()+arb(3).sqrt()*err
            keep=ds < arb(2)/3
            if keep:worst=max(worst,ds.upper())
        if not keep:
            assert d>=2
            half=d//2
            stack.extend((i+a,j+b,half) for a in (0,half) for b in (0,half))
    print("实际零点、覆盖方格数、行列式线对共用参考的主角正弦上界",
          n,visited,worst.str(12),flush=True)


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=128
    K,W,G=ingredients(shell,22)
    S=s_matrix(G,arb(1)/9,arb(3)/11)
    anchor,_=top_frame(K-W*S*adj(W))
    for n in (21,22):
        certify(shell,n,anchor)
    phi=(arb(4)/5).asin()
    assert 3*phi < arb.pi()
    print("每条行列式重叠边相角模的统一上界",phi.str(12))


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
