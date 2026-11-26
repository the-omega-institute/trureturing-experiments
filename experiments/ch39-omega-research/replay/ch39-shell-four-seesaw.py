"""Certified upper bounds on the seesaw gap for every positive mixture of the four shells."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import arb, arb_mat, acb

m=run_path(str(Path(__file__).with_name("ch39-shell-third-shell.py")))
dot,cross=m["dot"],m["cross"]
orbits={17:((0,1,4),(2,2,3)),867:((1,5,29),(13,13,23)),
        26:((0,1,5),(1,3,4)),29:((0,2,5),(2,3,4))}
eye=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])


def projection(v):
    return arb_mat([[arb(x*y)/dot(v,v) for y in v] for x in v])


def matrices(n):
    h=acb.zeta_zero(n).imag
    result=[]
    for r,(v,w) in orbits.items():
        points=sorted(set(permutations(w)),key=lambda p:-dot(v,p))
        mat=projection(v)
        for p,f in zip(points[1:3],m["shells"][r][2]):
            mat+=(2*h*f).sin()**2*projection(cross(cross(v,p),v))
        result.append(mat)
    return result


def array(M):
    return np.array([[float(M[i,j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


def positive(M):
    if M[0,0]<=0:
        return False
    return (M[0,0]>0 and M[0,0]*M[1,1]-M[0,1]*M[1,0]>0
            and (M.nrows()==2 or M.det()>0))


def spectra_at_vertex(mat):
    approx=array(mat)
    # Float values propose an upper bound only. Sylvester tests prove it.
    cap=arb(float(np.linalg.eigvalsh(approx)[-1])+1.e-8)
    assert positive(eye*cap-mat)
    return mat,approx,cap


def certify(n,limit):
    original=matrices(n)
    # Each simplex vertex stores its affine mass matrix.
    simplex=[spectra_at_vertex(a) for a in original]
    leaves=[]
    def visit(vertices,depth):
        assert depth<40
        mats,arrays,ceilings=zip(*vertices)
        center=sum(arrays)/4
        vec=np.linalg.eigh(center)[1]
        # These integer vectors choose one Rayleigh direction and a two-dimensional test space.
        vec=np.rint(vec*10**7).astype(int)
        q=arb_mat([[int(x)] for x in vec[:,0]])
        V=arb_mat([[int(x) for x in row] for row in vec[:,1:]])
        gram=V.transpose()*V
        assert positive(gram)
        A=max(((q.transpose()*M*q)[0,0]/(q.transpose()*q)[0,0]).upper() for M in mats)
        W=np.linalg.qr(vec[:,1:].astype(float),mode="reduced")[0]
        B=arb(float(min(np.linalg.eigvalsh(W.T@P@W)[0] for P in arrays))-1.e-8)
        assert all(positive(V.transpose()*M*V-gram*B) for M in mats)
        C=max(ceilings)
        if B>0 and C>A and A<B:
            bound=(A/B)**2*(C*C-B*B)/(C*C-A*A)
            if bound<limit:
                leaves.append(bound.upper())
                return
        # Bisect the longest edge in the affine image; both child simplices cover the parent.
        i,j=max(combinations(range(4),2),key=lambda ij:np.linalg.norm(arrays[ij[0]]-arrays[ij[1]]))
        mid=spectra_at_vertex((mats[i]+mats[j])/2)
        for index in (i,j):
            new=list(vertices)
            new[index]=mid
            visit(new,depth+1)
    visit(simplex,0)
    print("zero",n,"all nonnegative shell weights, upper bound",max(leaves).str(13),"simplices",len(leaves))
    return original


def main():
    for n,cap in ((25,".019"),(255,".024"),(719,".013")):
        certify(n,arb(cap))


if __name__=="__main__":
    main()
