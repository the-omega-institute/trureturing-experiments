"""Rotation selected by matching the two nine-point configurations and the cost of CP restoration."""
from itertools import product,permutations
from math import prod
from fractions import Fraction as F
from flint import arb,fmpz_mat,arb_mat

digits=(1,5,13,23,29)
def dots(x,y):
    return sum(a*b for a,b in zip(x,y))
def cross(x,y):
    return [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]
def unit(v):
    return [arb(x)/arb(dots(v,v)).sqrt() for x in v]
def shell(r):
    from math import isqrt
    return sorted((v for v in product(range(isqrt(r)+1),repeat=3) if dots(v,v)==r),
                  key=lambda v:(0 not in v,v))
sites= shell(17)
mapped=[tuple(digits[j] for j in v) for v in sites]
den=arb(17*867).sqrt()
def moment(ids):
    return fmpz_mat([[sum(sites[k][i]*mapped[k][j] for k in ids) for j in range(3)]
                    for i in range(3)])
C=moment(range(9))
transverse=C[0,0]-C[0,1]
longitudinal=C[0,0]+2*C[0,1]
assert all(C[i,j]==(C[0,0] if i==j else C[0,1]) for i,j in product(range(3),repeat=2))
for idx in (range(6),range(6,9)):
    mat=moment(idx)
    assert mat[0,0]>mat[0,1]>0
ug,vg=unit(sites[0]),unit(mapped[0])
cs=dots(ug,vg)
angle=cs.acos()
# CP is restored for generic nonzero amplitudes exactly when the two axial planes coincide.
# Minimize the labeled-shell least-squares cost under R*ug = +/- vg by reducing to a constrained O(2).
def plane_basis(n):
    u=unit(cross(n,[1,0,0]))
    return arb_mat([u,cross(n,u),n]).transpose()
A,B=plane_basis(ug),plane_basis(vg)
D=B.transpose()*arb_mat(C.tolist())*A
vals=[]
for sign in (1,-1):
    # orientation of the 2x2 block equals sign so that the full matrix is proper.
    a=D[0,0]+sign*D[1,1]
    b=D[1,0]-sign*D[0,1]
    vals.append(2*(arb(int(2*transverse+longitudinal))-sign*D[2,2]-(a*a+b*b).sqrt())/den)
assert 0<vals[0]<vals[1]
print("raw transverse and diagonal moments",transverse,longitudinal)
print("angle between shell basepoint radii",angle.str(13),
      "CP-restoring least extra sum of squared differences",vals[0].str(13))
for r,minimum in ((17,0),(867,145)):
    points=shell(r)
    winners=[tuple(sorted(v)) for v in points if prod(v)==minimum]
    assert len(winners)==6 and len(set(winners))==1 and min(map(prod,points))==minimum
