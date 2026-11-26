"""Family gauge thresholds if the two aligned spin-one frames provide their own inertia."""
from pathlib import Path
from fractions import Fraction as F
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-alignment-quartic.py")))
q,pi=m["q"],m["pi"]
f=m["fmin"].lower()
gu=(4*pi/q("40.9906")).sqrt()
maxN=515//9
# Defining-representation index of the family SU(3) subgroup (fundamental has level one).
reps=[]
for p in range(maxN):
    for r in range(p+1):
        dim=(p+1)*(r+1)*(p+r+2)//2
        if dim>maxN:
            break
        index=F(dim*(p*p+p*r+r*r+3*p+3*r),12)
        assert index.denominator==1
        if index:
            reps.append((dim,int(index)))
best=[0]*(maxN+1)
for n in range(1,maxN+1):
    best[n]=max([best[n-1]]+[idx+best[n-dim] for dim,idx in reps if dim<=n])
assert maxN==57 and best[maxN]==495
lower=f*gu*(q(9)/best[maxN]).sqrt()
print("family coupling upper to keep every vector below Q",(1/f).upper().str(12))
print("required family/color level ratio",(gu**2*f**2).lower().str(12))
print("mass/Q lower for at least three family vectors with gF=gU, or with allowed SU(N) family level",
      (f*gu).lower().str(12),lower.lower().str(12))
