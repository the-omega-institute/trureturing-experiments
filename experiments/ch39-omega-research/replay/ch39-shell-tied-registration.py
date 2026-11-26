"""Unlabeled registration of the two twelve-site shells, including orbit-weighted stiffness."""
from fractions import Fraction
from itertools import permutations
from math import isqrt
from flint import fmpz_mat

outer=(0,1,5)
inner=(1,3,4)
radii=26,29
# Start with non-neighboring measurements on the widest orbit.
order=((0,1,5),(5,1,0),(1,5,0),(0,5,1),(5,0,1),(1,0,5))
source=order+tuple(sorted(set(permutations(inner))))
target=tuple(sorted(set(permutations((0,2,5)))|set(permutations((2,3,4)))))
id4=fmpz_mat([[int(i==j) for i in range(4)] for j in range(4)])


def K(x,y):
    C=[[a*b for b in y] for a in x]
    tr=sum(C[i][i] for i in range(3))
    v=[C[(i+1)%3][(i+2)%3]-C[(i+2)%3][(i+1)%3] for i in range(3)]
    return fmpz_mat([[tr,*v]]+[[v[i],*[C[i][j]+C[j][i]-int(i==j)*tr for j in range(3)]]
                                for i in range(3)])


def pd(M):
    a=M[0,0]
    if a<=0:
        return False
    b=a*M[1,1]-M[0,1]**2
    if b<=0:
        return False
    c=b*M[2,2]-a*M[1,2]**2-M[1,1]*M[0,2]**2+2*M[0,1]*M[1,2]*M[0,2]
    return c>0 and M.det()>0


def diag_map(x):
    i=0 if 0 in x else 1
    shape=(outer,inner)[i]
    dest=((0,2,5),(2,3,4))[i]
    return tuple(dest[shape.index(xi)] for xi in x)


def verify(outer_weight,inner_weight):
    weights=6*[outer_weight]+6*[inner_weight]
    desired=sum(w*sum(a*b for a,b in zip(x,diag_map(x))) for x,w in zip(source,weights))
    # Resolve assignments whose correlation is within one integer unit of this witness.
    cut=desired-1
    blocks=[[w*K(x,y) for y in target] for w,x in zip(weights,source)]
    tail=[sum(weights[i:]) for i in range(13)]
    bound=[isqrt((w*1000)**2*26*29)+1 for w in tail]
    full=set(range(12))
    final=[]
    nodes=0
    # Cyclic coordinate permutations on the target are proper rotations.
    first={j for j,y in enumerate(target) if y==min(y,y[1:]+y[:1],y[2:]+y[:2])}
    def visit(seq,S):
        nonlocal nodes
        nodes+=1
        k=len(seq)
        if pd((1000*cut-bound[k])*id4-1000*S):
            return
        if k==12:
            if not pd(cut*id4-S):
                final.append((seq,S))
            return
        for j in (first if not k else full.difference(seq)):
            visit(seq+(j,),S+blocks[k][j])
    visit((),fmpz_mat(4,4))
    expected=tuple(target.index(diag_map(x)) for x in source)
    assert expected[0] in first
    opt=[]
    for seq,S in final:
        if seq==expected:
            assert S[0,0]==desired and all(S[0,i]==0 for i in (1,2,3))
            assert pd(desired*fmpz_mat([[int(i==j) for j in range(3)] for i in range(3)])
                      -fmpz_mat([[S[i,j] for j in (1,2,3)] for i in (1,2,3)]))
            opt.append(seq)
        else:
            assert pd(desired*id4-S)
    assert opt==[expected]
    s=39*outer_weight+9*inner_weight
    ell=84*outer_weight+144*inner_weight
    assert desired==2*s+ell
    print("orbit weights",outer_weight,inner_weight,"maximum trace",desired,
          "squared frequency ratio",Fraction(2*s,s+ell),"visited",nodes,"above cutoff",len(final))
    return Fraction(2*s,s+ell)


if __name__=="__main__":
    assert verify(1,1)==Fraction(8,23)
    assert verify(15,7)==Fraction(4,9)
