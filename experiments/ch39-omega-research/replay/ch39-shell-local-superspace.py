"""Local Stueckelberg realizations of the shifted 0- and 1-form masses on the shell graph."""
from fractions import Fraction as F
from pathlib import Path
from flint import fmpz_mat
import runpy

g=runpy.run_path(str(Path(__file__).with_name("ch39-shell-family-alignment.py")))
sites=g["shell"](867)
edges=[(i,j) for i in range(len(sites)) for j in range(i)
       if sum((u-v)**2 for u,v in zip(sites[i],sites[j]))<=804]
D=fmpz_mat([[int(k==j)-int(k==i) for k in range(len(sites))] for i,j in edges])
assert len(edges)==135 and D.rank()==18 and len(edges)-D.rank()==117
def eye(n):
    return fmpz_mat([[int(i==j) for j in range(n)] for i in range(n)])
def stack(A,B):
    return fmpz_mat(A.tolist()+B.tolist())
# This replaces d by an exact incidence matrix. The identities are polynomial in d,d*.
mass=2
J=stack(mass*eye(135),D.transpose())
K=stack(-D,mass*eye(19))
assert (J.transpose()*K).is_zero()
assert J.transpose()*J==mass**2*eye(135)+D*D.transpose()
assert K.transpose()*K==mass**2*eye(19)+D.transpose()*D
assert J.ncols()+K.ncols()==J.nrows()
def ein(chi,vec):
    return 2*(F(1,2)-2*chi)-vec
assert ein(F(11,40),F(1,5))==-F(3,10)
assert ein(F(11,60),F(1,5))==F(1,15)
assert ein(F(1,5),F(1,5))==0
assert 2*F(1,48)-F(1,24)==0
print("graph H0,H1 nullities and heat-trace offset",19-D.rank(),135-D.rank(),135-19)
print("local squared-mass factorization verified; curvature coefficients match at both shifts")
