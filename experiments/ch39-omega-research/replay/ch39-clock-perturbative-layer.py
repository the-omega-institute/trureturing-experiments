"""A uniform bound for layers whose weak-doublet lift meets the matching-scale charge criterion."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import runpy
from flint import arb

study=runpy.run_path(str(Path(__file__).with_name("ch39-clock-variable-doublet-budget.py")))
m,ball=study["m"],study["ball"]
def charges(t):
    p,q,u,v,r,s=t
    return (F(p*p+q*q+p*q+3*p+3*q,3),F(3*u*u+2*u*v+3*v*v+12*u+12*v,8),
            F(3,4),F(r*(r+2),4),F(s*(s+2),4))
def min_charge(t):
    x,c,l,r,s=charges(t)
    return ball(c+l)+(ball(F(16,3)*x).sqrt()+ball(r).sqrt()+ball(F(4,9)*s).sqrt())**2
maximum=arb("40.9906")
limits=[]
for i in range(6):
    t=[0]*6
    while min_charge(t)<maximum:
        t[i]+=1
    assert min_charge(t)>maximum
    limits.append(t[i])
# The lower bounds on the Casimirs increase with each Dynkin label. Prune complete boxes.
weights=[F(273,433),0,0,F(1456,4763),F(304,4763)]
opt=F(10920,4763)
assert sum(weights)==1
costs=[]
def visit(prefix):
    t=prefix+(0,)*(6-len(prefix))
    lower=min_charge(t)
    if lower>maximum:
        return
    assert lower<maximum
    if len(prefix)<6:
        for i in range(limits[len(prefix)]):
            visit(prefix+(i,))
    elif m["target"](*t):
        v=study["coeff"](*t)
        costs.append(sum(w*b for w,b in zip(weights,v)))
visit(())
assert costs and min(costs)==opt
c=[6*w/F(11*n) for w,n in zip(weights,m["ca"])]
budget=2*arb.pi()*maximum*ball(c[1]+max(F(3,16)*c[0],c[3],F(9,4)*c[4]))
def integrated(x):
    return ball(opt-2)*study["kappa"]*(x-1)-ball(270*opt+78)*x.log()
lo,hi=arb(1),arb(559)
for _ in range(40):
    mid=(lo+hi)/2
    if integrated(mid)>budget:
        hi=mid
    else:
        lo=mid
assert integrated(hi)>budget
print("exclusive label ranges",limits,"irreps satisfying relaxed charge and target",len(costs))
print("dual cost per doublet",opt,"two-tower pole upper E/Q",hi.upper().str(13))
