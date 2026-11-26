"""Global smooth square-root sizes for complex and quaternionic shell holonomies, including -1."""
from pathlib import Path
from runpy import run_path

from flint import fmpq as Q, fmpq_mat as Mat
from sympy import Matrix as SymMatrix

orient=run_path(str(Path(__file__).with_name("ch39-shell-nine-phase-probes.py")))["orientation"]


def eye(n):
    a=Mat(n,n)
    for i in range(n):
        a[i,i]=1
    return a


def left(q):
    if len(q)==2:
        a,b=q
        return Mat([[a,-b],[b,a]])
    a,b,c,d=q
    return Mat([[a,-b,-c,-d],[b,a,-d,c],[c,d,a,-b],[d,-c,b,a]])


def main():
    edges=orient()
    ordinary=SymMatrix.zeros(48,24)
    for i,(a,b) in enumerate(edges):
        ordinary[i,a],ordinary[i,b]=-1,1
    tree=list(ordinary.T.rref()[1])
    chords=[e for e in range(48) if e not in tree]
    for dim in (2,4):
        unit=[1]+[0]*(dim-1)
        hol={e:unit[:] for e in range(48)}
        hol[chords[0]]=[-1]+[0]*(dim-1)
        hol[chords[1]]=[Q(3,5),Q(4,5)]+[0]*(dim-2)
        if dim==4:
            hol[chords[2]]=[Q(5,13),0,Q(12,13),0]
            assert left(hol[chords[1]])*left(hol[chords[2]])!=left(hol[chords[2]])*left(hol[chords[1]])
        D=Mat(dim*48,dim*24)
        T=Mat(dim*23,dim*23)
        C=Mat(dim*25,dim*23)
        b=Mat(dim*25,1)
        for seq,e in enumerate(tree+chords):
            v,w=edges[e]
            for j in range(dim):
                for k in range(dim):
                    for site,entry in ((v,-int(j==k)),(w,left(hol[e])[j,k])):
                        D[dim*seq+j,dim*site+k]=entry
                        if site:
                            if seq<23:
                                T[dim*seq+j,dim*(site-1)+k]=entry
                            else:
                                C[dim*(seq-23)+j,dim*(site-1)+k]=entry
            if seq>=23:
                for j in range(dim):
                    b[dim*(seq-23)+j,0]=hol[e][j]-unit[j]
        A=C*T.inv()
        H=(eye(dim*25)+A*A.transpose()).inv()
        s=(b.transpose()*H*b)[0,0]
        # Fix vertex 0 to the scalar one and solve the grounded least squares.
        grd=Mat([[D[i,j] for j in range(dim,24*dim)] for i in range(48*dim)])
        root=Mat([[D[i,0]] for i in range(48*dim)])
        raw=root-grd*(grd.transpose()*grd).solve(grd.transpose()*root)
        assert (raw.transpose()*raw)[0,0]==s and s>0
        dim_im=25*(dim-1)
        features=Mat(dim_im+1,1)
        combine=Mat(dim*25,dim_im+1)
        relation=Mat(1,dim_im+1)
        deficit=sum(1-hol[e][0] for e in chords)
        features[dim_im,0]=deficit
        relation[0,dim_im]=1
        for pos,e in enumerate(chords):
            d=1-hol[e][0]
            den=2+deficit-d
            beta=Q(d)/den
            relation[0,dim_im]-=beta
            combine[dim*pos,dim_im]=-beta
            for j in range(dim-1):
                v=hol[e][j+1]
                col=(dim-1)*pos+j
                alpha=Q(v)/den
                combine[dim*pos,col]=-alpha
                combine[dim*pos+j+1,col]=1
                features[col,0]=v
                relation[0,col]=-alpha
        assert combine*features==b and (relation*features)[0,0]==0
        metric=combine.transpose()*H*combine+relation.transpose()*relation
        assert metric.det()>0 and (features.transpose()*metric*features)[0,0]==s
        rows_free=23+(dim_im+1)//dim
        assert dim*(rows_free-24)<dim_im<=dim*(rows_free-23)
        print("Real dimension of each coefficient",dim,"flux directions",dim_im,
              "global analytic rows",rows_free,"example exact grounded residual",s)


if __name__=="__main__":
    main()
