"""Screening lower bound for either sign of a color- and weak-singlet unit Weyl charge."""
from fractions import Fraction as F
from pathlib import Path
from flint import arb
import runpy

study=runpy.run_path(str(Path(__file__).with_name("ch39-clock-variable-doublet-budget.py")))
m,ball=study["m"],study["ball"]
def dim(t):
    p,q,u,v,r,s=t
    return m["xdim"](p,q)*m["cdim"](u,v)*(r+1)*(s+1)
def weyl(t):
    relaxed=study["weighted"](study["coeff"](*t,relaxed=True))
    return 2*relaxed  # half the weak dimension, Weyl factor (2/3)/(1/6)
def total_unit(t):
    p,q,u,v,r,s=t
    return sum(i-j+1 for i in range(q,p+q+1) for j in range(q+1)
               for z in range(-r,r+1,2) for a in range(-s,s+1,2)
               if abs(4*(3*(i+j)-2*(p+2*q))+3*(v-u)+3*z+2*a)==6)
def lower(t):
    # One color singlet in each admissible SU(4) irrep; divide by its maximum possible multiplicity.
    return weyl(t)*m["cdim"](t[2],t[3])/dim(t)
origin=(0,0,0,0,2,0)
best=weyl(origin)/total_unit(origin)
winners=[]
def search(prefix):
    global best,winners
    if len(prefix)<6:
        n=0
        while lower(prefix+(n,)+(0,)*(5-len(prefix)))<=best:
            search(prefix+(n,))
            n+=1
        return
    num=total_unit(prefix)
    if num and weyl(prefix)<=best*num:
        value=weyl(prefix)/num
        if value<best:
            best,winners=value,[]
        winners.append(prefix)
search(())
print("per charged left Weyl component",best,float(best),"winners",winners)
assert best>0
for y,ns,extra in ((4,6,F(39585,119207)),(5,4,F(10920,10837))):
    # The right and left graphs have identical optical length but different thresholds.
    k=study["opt"]+extra+ns*best
    def expenditure(x):
        return ball(k-2)*study["kappa"]*(x-1)-ball(270*k+78)*x.log()
    print("total minimum against the chosen dual weight",y,k,float(k))
    if k<=2:
        continue
    lo,hi=arb(1),arb(559)
    assert expenditure(hi)>study["budget"]
    for _ in range(40):
        mid=(lo+hi)/2
        if expenditure(mid)>study["budget"]:
            hi=mid
        else:
            lo=mid
    print("two gauge towers per factor: pole E/Q upper",hi.upper().str(13))
