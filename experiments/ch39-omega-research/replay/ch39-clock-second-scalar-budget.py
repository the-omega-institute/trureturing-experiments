"""Common separating inequality for both large scalars of each right shell layer."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

study=runpy.run_path(str(Path(__file__).with_name('ch39-clock-variable-doublet-budget.py')))
weights,base,ball,kappa,budget=[study[key] for key in ('weights','opt','ball','kappa','budget')]
def relaxed_cost(labels,real):
    p,q,u,v,r,s=labels
    # weak singlet halves the dimension, and its SU(2)_L component has zero dual weight.
    cost=study["weighted"](study["coeff"](*labels,relaxed=True))/2
    return cost if real else 2*cost

def singlet(y):
    best=relaxed_cost((0,0,0,0,2*y,0),True)
    winners=[]
    def search(prefix):
        nonlocal best,winners
        t=prefix+(0,)*(6-len(prefix))
        floor=relaxed_cost(t,True)
        if floor>best:
            return
        if len(prefix)<6:
            n=0
            while relaxed_cost(prefix+(n,)+(0,)*(5-len(prefix)),True)<=best:
                search(prefix+(n,))
                n+=1
            return
        p,q,u,v,r,s=t
        real=p==q and u==v and (r+s)%2==0
        cost=relaxed_cost(t,real)
        if cost>best:
            return
        hs={3*(i+j)-2*(p+2*q) for i in range(q,p+q+1) for j in range(q+1)}
        if any(abs(4*h+3*(v-u)+3*z+2*w)==6*y for h in hs
               for z in range(-r,r+1,2) for w in range(-s,s+1,2)):
            if cost<best:
                best,winners=cost,[]
            winners.append(t)
    search(())
    return best,winners
for y in (4,5):
    b,winners=singlet(y)
    k=base+b
    print("Y",y,"singlet dual minimum",b,winners,"combined",k,float(k))
    for copies in (1,2):
        if k<=copies:
            continue
        def expenditure(x):
            return ball(k-copies)*kappa*(x-1)-ball(270*k+38*copies+2)*x.log()
        lo,hi=arb(1),arb(559)
        assert expenditure(hi)>budget
        for _ in range(40):
            mid=(lo+hi)/2
            if expenditure(mid)>budget:
                hi=mid
            else:
                lo=mid
        print("gauge tower copies",copies,"max E/Q",hi.upper().str(12))
