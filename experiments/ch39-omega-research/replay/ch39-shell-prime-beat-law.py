"""Character lattice of the finite-window prime terms and the leading covariance law."""
from flint import fmpz_mat,arb
from pathlib import Path
import runpy
from itertools import combinations,product

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-prime-covariance.py")))
local=m["m"]
l,r=m["l"],m["r"]
directions=[(1,0),(0,1),(-1,1),(1,1)]
rows=[]
gaps=[]
def record(left,right,sign,prime):
    value=sum(i*f for i,f in zip(left+right,[l[0],l[1],r[0],r[1]]))
    assert sign*value>0
    gap=sign*value-arb(prime).log()
    assert abs(gap)<local["B"] and not gap.contains(0)
    rows.append([sign*c for c in left+right]+[-int(prime==p) for p in (2,3,5,19)])
    gaps.append(gap)
record((0,1),(0,0),1,2)
for i,j,k in ((1,1,1),(1,3,-1),(3,3,1),(3,2,-1),(2,0,-1),(2,1,1),(2,3,-1)):
    sign=1 if l[i]+k*r[j]>0 else -1
    record(directions[i],tuple(k*b for b in directions[j]),sign,m["hits"][i,j,k])
M=fmpz_mat(rows)
assert M.rank()==6
relations=fmpz_mat([[-1,0,-1,0,0,1,1,0],[1,-1,0,0,0,-1,0,1]])
assert (relations*M).is_zero()
assert relations[0,6]==1 and relations[1,7]==1 and relations[0,7]==relations[1,6]==0
assert fmpz_mat(rows[:6]).rank()==6
# Leading coefficient of log(T/(2pi))*Cov(t17,t867), ignoring the sign of the cosine.
gap=abs(gaps[1])
amp=arb(5).log()/(8*arb(5).sqrt())*local["density"](gap,gap<2*local["a"])/local["f0"]
assert arb(".08")<amp<arb(".10")
print("eight signed detunings",*[x.str(10) for x in gaps])
print("independent detunings",M.rank(),"scaled covariance amplitude",amp.str(15),
      "limiting standard deviation",(amp/arb(2).sqrt()).str(14))
