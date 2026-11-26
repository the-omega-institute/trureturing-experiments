"""Sharp kinetic ratio for three net copies in every color/lepton direction, including face ties."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
from runpy import run_path
from flint import arb

home=Path(__file__).resolve().parent
d=run_path(str(home/"ch39-shell-octahedral-determinant.py"))


def costs(a,b):
    allrows={k:(d["traceless_cost"](4,k-a),d["traceless_cost"](3,-k-b),F(k*k))
             for k in range(-6,7)}
    return {k:v for k,v in allrows.items() if not any(
        v!=w and all(x>=y for x,y in zip(v,w)) for w in allrows.values())}


def color(k,base):
    # A shortest log contains two adjacent lifts; the first three entries must be equal.
    ans=[p for p in range(-2,3) if abs(k-base-4*p)<=1]
    return None if not ans else F(base,4)+ans[0]


def extreme(ineq):
    def cross(v,w):
        return tuple(v[(i+1)%3]*w[(i+2)%3]-v[(i+2)%3]*w[(i+1)%3] for i in range(3))
    rays=set()
    for a,b in combinations(ineq,2):
        p=cross(a,b)
        if not sum(p):
            continue
        p=tuple(j/sum(p) for j in p)
        if all(sum(x*y for x,y in zip(p,row))>=0 for row in ineq):
            rays.add(p)
    if not rays or any(all(p[j]==0 for p in rays) for j in range(3)):
        return set()
    return rays


def surface():
    caps,sides=costs(1,0),costs(-1,-1)
    lower=F(1000)
    sharp=[]
    for c,s in product(combinations_with_replacement(caps,2),combinations_with_replacement(sides,6)):
        K=sum(c+s)
        colors=[color(k,base) for group,base in ((c,1),(s,-1)) for k in group]
        if None in colors or abs(K)!=12 or 4*sum(colors)!=K:
            continue
        inequalities={(F(1),0,0),(0,F(1),0),(0,0,F(1))}
        for table,chosen in ((caps,c),(sides,s)):
            for k in set(chosen):
                for v in table.values():
                    inequalities.add(tuple(x-y for x,y in zip(v,table[k])))
        rays=extreme(inequalities)
        if rays:
            lo=min(p[0]/p[1] for p in rays if p[1])
            lower=min(lower,lo)
            for p in rays:
                if p[1] and p[0]/p[1]==F(2,3) and p[2]:
                    sharp.append((c,s,p[2]/p[1]))
    assert lower==F(2,3) and sharp
    print("sharp necessary g_F^2/g_4^2",lower,"first tie witness",sharp[0],"patterns at boundary",len(sharp))
    return lower


if __name__=="__main__":
    ratio=surface()
    budget=run_path(str(home/"ch39-shell-point-family-budget.py"))
    gu,gmax=budget["gu"],budget["gmax"]
    g=gu*budget["q"](ratio).sqrt()
    print("at least two family vector masses /Q above",(g/gmax).str(13))
    print("induced inverse Newton fraction if vectors <=Q, upper",budget["upper"](budget["cap"]/g**2).upper().str(13))
