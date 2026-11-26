"""Neutral-companion contribution to the finite Einstein term on the joint solution."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-neutral-clock-consecutive.py"))
q,r=m["q"],m["r"]
R,_=m["rotations"](719)
lo=q(".25")
step=q(2)**(q(1)/16)
# At positive shift K0 >= L + z I. The ordinary scalar block has gap >8 above
# the one constant mode. The twisted endpoint norm grows at most (1+beta)z.
src=m["pot"]["short"]
cap=src["base"]["arc_budget"]()[1].upper()
beta=9*q("1.018")**2*m["sh"]["rho"]*m["sh"]["theta"](804)
top=cap+(1+beta)*lo
upper=lo*(57*top.log()-54*q(8).log()-3*lo.log()+3)
lower=q(0)
for _ in range(16*32):
    hi=lo*step
    lower+=(hi-lo)*m["point"](hi,False,R).lower()
    upper+=(hi-lo)*m["point"](lo,True,R).upper()
    lo=hi
degree=cap/2
upper+=608*degree**3/lo**2
value=lower.lower().union(upper.upper())/r**2
assert value>0
inc=value/(4*arb.pi())
ql=1000*q("26.8220").exp()
qh=1000*q("27.0455").exp()
Mp=q("1.220890e19")
ratio=((2*q(".05")/(3*arb.pi())+inc.lower())*(ql/Mp)**2).lower().union(
       ((2*q("1.99779")/(3*arb.pi())+inc.upper())*(qh/Mp)**2).upper())
print("neutral dimensionless Einstein heat moment",value.lower().str(12),value.upper().str(12))
print("total induced fraction of measured inverse G",ratio.lower().str(12),ratio.upper().str(12))
assert ratio<q("1e-9")
