"""Allow a finite shift of the hypercharge matching while retaining exact color-weak matching."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
from collections import Counter
from flint import fmpz_mat
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-simple-perturbative.py")))
ceiling=m["ceiling"]
tolerance=F(1,20)  # absolute shift of g_U^2/g_Y^2, i.e. 3% of its canonical value
systems=[]
for c in range(1,ceiling//11+1):
    upper=F(c)*(1+tolerance)
    for x in range(1,int(3*upper/16)+1):
        for r in range(1,int(upper-16*F(x,3))+1):
            for s in range(1,int((upper-16*F(x,3)-r)*9/4)+1):
                shift=F(48*x+9*r+4*s-9*c,9*c)
                if abs(shift)<=tolerance and shift:
                    systems.append((x,c,c,r,s))
dimension=max(ceiling//s[1] for s in systems)
context=m["blocks"].__globals__
context["max_dim"]=dimension
context["reps"]=tuple(m["simple"](n) for n in (3,4,2,2,2))
samples=m["samples"]
y=[m["light"](t) for t in samples]
checked=Counter()
survive=[]
for target in systems:
    x,c,_,r,s=target
    nmax=ceiling//c
    B=m["blocks"](target,nmax)
    count,multisets=m["suffix_solver"](B,target,nmax)
    if not count(target,nmax,0):
        continue
    tr2=[c*(a*a+b*b+(a+b)**2+2*w*w)+(96*x+12*c+18*r+8*s)*h*h for a,b,w,h in samples]
    for prefix in multisets(target,nmax,0,()):
        for orient in product(*(m["orientations"](B[k][2]) for k in prefix)):
            rows=[m["moment"](blk) for blk in orient]
            matrix=fmpz_mat([[sum(z[i][1] for z in rows),tr2[i]*sum(z[i][0] for z in rows),y[i]]
                             for i in range(len(samples))])
            checked[target]+=1
            if matrix.rank()<3 and matrix.rank()==fmpz_mat([[matrix[j,0],matrix[j,1]] for j in range(len(samples))]).rank():
                survive.append((target,orient))
    count.cache_clear()
print("distinct noncanonical level systems with a defining representation",len(checked))
print("checked oriented decompositions away from exact matching",sum(checked.values()),
      "passing seven samples",len(survive))
assert len(checked)==61 and sum(checked.values())==174814 and not survive
