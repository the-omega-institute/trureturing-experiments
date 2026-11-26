"""第一与第二枚实际零点的矩阵弦，在指定两条边全非负权象限上的三四带隙。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


EDGES=((0,2),(5,14))
UNIT=65536


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    sh=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=112
    heights=[acb.zeta_zero(i).imag for i in (1,2)]
    for q in (1,-4,2,-3,6):
        A,Z=[sh["laplacian"](q,h) for h in heights]
        D=Z-A
        norm=max(sum(abs(D[i,j]) for j in range(19)).upper() for i in range(19))
        columns=[]
        phases=[]
        for i,j in EDGES:
            v=acb_mat(19,1)
            v[i,0]=1
            v[j,0]=A[j,i]
            assert abs(A[j,i]).contains(1) and A[j,i].overlaps(Z[j,i])
            phases.append(-A[i,j])
            columns.append(v*adj(v))
        @cache
        def compression(infinite):
            used=set()
            vectors=[]
            for k in infinite:
                i,j=EDGES[k]
                v=[acb(0)]*19
                v[i],v[j]=phases[k]/arb(2).sqrt(),1/arb(2).sqrt()
                used.update((i,j))
                vectors.append(v)
            for i in range(19):
                if i not in used:
                    v=[acb(0)]*19
                    v[i]=1
                    vectors.append(v)
            return acb_mat(vectors).transpose()
        @cache
        def corner(t,x,y):
            h=A+(arb(t)/(2*UNIT))*D
            infinite=[]
            for k,p in enumerate((x,y)):
                if p==UNIT:
                    infinite.append(k)
                else:
                    h+=columns[k]*(arb(8*p)/(UNIT-p)-1)
            if infinite:
                basis=compression(tuple(infinite))
                h=adj(basis)*h*basis
            ev=sorted((z.real for z in h.eig(multiple=True)),key=float)
            assert ev[1]<ev[2]<ev[3]<ev[4]
            return ev[2],ev[3]
        # 时间中点编码为上下两个整数格号之和。
        stack=[(0,UNIT,0,UNIT,0,UNIT)]
        checked=0
        minimum=arb(100)
        while stack:
            t,T,x,X,y,Y=stack.pop()
            gap=corner(t+T,x,y)[1]-corner(t+T,X,Y)[0]-norm*(T-t)/UNIT
            if gap>arb(".0001"):
                checked+=1
                minimum=min(minimum,gap.lower())
                continue
            widths=[T-t,X-x,Y-y]
            axis=max(range(3),key=lambda i:widths[i])
            assert widths[axis]>1,(q,t,T,x,X,y,Y,gap)
            box=[t,T,x,X,y,Y]
            a,b=2*axis,2*axis+1
            mid=(box[a]+box[b])//2
            c=box.copy()
            c[a]=mid
            box[b]=mid
            stack.extend((tuple(c),tuple(box)))
        print("荷、整弦全非负双权认证盒数、统一谱隙严格下界",q,checked,minimum.str(12),flush=True)
        corner.cache_clear()


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
