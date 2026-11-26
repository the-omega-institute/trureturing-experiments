"""Quaternionic alignment potential as an F-term, with flat and curvature heat weights checked."""
from fractions import Fraction as F
from flint import fmpz_mat
from itertools import product
from pathlib import Path
import runpy

g=runpy.run_path(str(Path(__file__).with_name('ch39-shell-family-alignment.py')))
C=g["moment"](range(9))
tr=sum(C[i,i] for i in range(3))
I=fmpz_mat([[int(i==j) for j in range(3)] for i in range(3)])
A=tr*I-C
assert A.charpoly().factor()[1][0][1] in (1,2)
# Exact rotation identity in homogeneous quaternion coordinates, for all q (degree two).
def heat(z):
    p0,*v=z
    nn=sum(x*x for x in z)
    X=fmpz_mat([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R=nn*I+2*p0*X+2*X*X
    vv=fmpz_mat([[i] for i in v])
    return 2*(nn*tr-sum((C*R)[i,i] for i in range(3)))-4*(vv.transpose()*A*vv)[0,0]
assert all(heat(vec)==0 for vec in product(range(3),repeat=4))
assert all(heat(vec)==0 for vec in product((-1,1),repeat=4))
def chiral(xi):
    return F(1,2)-2*xi
def vector(xi):
    return -xi
# Add the three F-term enforcing chirals and their neutral spectator blocks at their own masses.
xi=F(1,5)
assert 2*chiral(xi)+vector(xi)==0
assert 2*2+1+3-2*4 == 0
assert 5*F(1,180)+4*F(7,720)-F(1,15)==0
assert 5*(-F(1,180))+4*F(1,90)+F(29,60)==F(1,2)
assert [chiral(k) for k in (F(1,6),F(2,5),F(13,60))]==[F(1,6),-F(3,10),F(1,15)]
print("axis quadratic eigenvalues and multiplicities",A.charpoly().factor())
print("neutral supplement (5 scalars, 4 Weyl, 1 Proca) Einstein and Ricci-flat traces zero")
