"""Minimum number of simple-source quadratic forms for the bare integral GS pair and the chiral metric."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import fmpz_mat

m=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-determinant.py")))
# Use coordinates (P,l,r), P=-ch2(A tensor B). Negative p1/2 pulls back to these
# classes for the spin lifts of both realified complex bundles.
X=(0,-1,1)
minus_Y=(1,1,1)
K1=(1,1,1)
K2=(1,0,2)
assert tuple(b-a for a,b in zip(K1,K2))==X and K1==minus_Y
polar=fmpz_mat([[X[i]*minus_Y[j]+minus_Y[i]*X[j] for j in range(3)] for i in range(3)])
assert polar.rank()==2
def face(k,kind,u,v):
    a,b=(1,0) if kind=="cap" else (-1,-1)
    return ((3*u+11*v)*m["traceless_cost"](4,k-a)+(4*u+11*v)*m["traceless_cost"](3,-k-b)
            +(u+5*v)*k*k/12)
def minima(kind,u,v):
    # The order-two determinant lift fixes the odd face sector.
    vals={k:face(k,kind,u,v) for k in (-5,-3,-1,1,3,5)}
    return [k for k,w in vals.items() if w==min(vals.values())]
for v,cap,side in ((0,[-3],[-5]),(F(1,4),[-3],[-5,-1]),(F(1,2),[-3],[-1])):
    assert minima("cap",F(1),v)==cap and minima("side",F(1),v)==side
for kind,chosen in (("cap",-3),("side",-1)):
    for k in (-5,-3,-1,1,3,5):
        if k!=chosen:
            assert face(k,kind,F(1),F(1,4))>=face(chosen,kind,F(1),F(1,4))
            assert face(k,kind,F(0),F(1))>face(chosen,kind,F(0),F(1))
assert face(-1,"side",F(1),F(0))-face(-5,"side",F(1),F(0))==F(2,3)
assert face(-1,"side",F(0),F(1))-face(-5,"side",F(0),F(1))==F(-8,3)
u1=u2=v=F(1,2)
color=3*(u1+u2)+11*v
left=u1+16*v
right=u1+2*u2+14*v
assert left==color==right==F(17,2)
print("rank of mixed polynomial",polar.rank(),"two sources select mixed index per color 9")
print("sharp third-source weight",F(1,4),"example local levels color=left=right",color,
      "third defining dimension",39+2*(16+14))
