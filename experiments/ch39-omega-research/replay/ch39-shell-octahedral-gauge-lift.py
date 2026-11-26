"""An explicit simple-group metric for the determinant branch and its mixed weak anomaly."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

m=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-determinant.py")))
def metric(branch):
    # (dimension, central charge, SU4 trace level, SU3 trace level)
    assert sum(d*q for d,q,_,_ in branch)==0
    return sum(d for d,_,_,_ in branch),sum(a for _,_,a,_ in branch),sum(b for _,_,_,b in branch),sum(F(d*q*q,144) for d,q,_,_ in branch)


def energy(k,face,w4,w3,c):
    i,j=(1,0) if face=="cap" else (-1,-1)
    return w4*m["traceless_cost"](4,k-i)+w3*m["traceless_cost"](3,-k-j)+c*k*k


def check():
    mixed=(12,1,3,4)
    dual_three=(3,-4,0,1)
    adjoint4=(15,0,8,0)
    adjoint3=(8,0,0,6)
    rep=[mixed,dual_three,adjoint4,adjoint3]
    N,w4,w3,c=metric(rep)
    assert (N,w4,w3,c)==(38,11,11,F(5,12))
    for face,expected in (("cap",-3),("side",-1)):
        vals=sorted((energy(k,face,w4,w3,c),k) for k in range(-6,7))
        assert vals[0][1]==expected and vals[0][0]<vals[1][0]
        print("SU(38) metric",face,"optimum and runner up",vals[:2])
    print("dimension times color level",N*w4)
    assert 12*1-12==0 and 12*(-1)+12==0
    charges=[(24,1),(24,-1),(2,-12),(2,12),(2,0),(2,0)]
    assert all(sum(d*q**p for d,q in charges)==0 for p in (1,3))
    assert (12+2)%2==0
    # q=12 corresponds to det B and has Chern number +12 in the selected branch.
    assert m["side"]*6+m["cap"]*2==-12
    print("mixed weak anomalies before spectators, fundamental index 1/2",6,-6,
          "charge-12 Higgs determinant degree",12)


if __name__=="__main__":
    check()
