"""Uniform singular-value separation for single-point, flavor-universal Yukawa evaluation."""
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
from runpy import run_path
from flint import arb

m=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-clock-modes.py")))
rows=range(3)
off=list(product(combinations(rows,2),repeat=2))
target=arb(".065")


def interval(cell):
    level,index=cell
    return arb(index).union(index+1)/2**level


@lru_cache(None)
def coeff(cell1,cell2):
    u=m["comps"](interval(cell1),m["ends"][0])
    v=m["comps"](interval(cell2),m["ends"][1])
    cs=[]
    for j in rows:
        a=[u[(j-k)%3]*m["zplus"]**k for k in rows]
        b=[v[(j-k)%3]*m["zminus"]**k for k in rows]
        c=[a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in rows]
        den=(abs(c[0])**2+abs(c[1])**2/2+abs(c[2])**2).sqrt()
        if not den>0:
            return None
        cs.append([v/den for v in c])
    return cs


@lru_cache(None)
def polys(cell1,cell2,reverse):
    cs=coeff(cell1,cell2)
    if cs is None:
        return None
    T=[arb(0) for _ in rows]
    terms=[]
    def entry(i,j):
        k=(j-i)%3
        return cs[j][k],2-k if reverse else k
    for i,j in product(rows,repeat=2):
        a,k=entry(i,j)
        T[k]+=abs(a)**2
    for (i,k),(j,l) in off:
        a,e=entry(i,j)
        b,f=entry(k,l)
        c,g=entry(i,l)
        d,h=entry(k,j)
        p,q=e+f,g+h
        if p==q:
            terms.append((abs(a*b-c*d),p,arb(0),q))
        else:
            terms.append((abs(a*b),p,abs(c*d),q))
    return T,terms


def prove(cells,reverse):
    s,t,r=cells
    data=polys(s,t,reverse)
    if data is not None:
        T,terms=data
        x=interval(r)
        powers=[arb(1)]+[x**k for k in range(1,5)]
        norm=sum(t*powers[2*k] for k,t in enumerate(T))
        minor=sum((abs(a*powers[p]-b*powers[q]).abs_lower())**2 for a,p,b,q in terms)
        if minor>3*target**2*norm**2:
            return 1
    i=min(range(3),key=lambda k:cells[k][0]+(1 if data is None and k==2 else 0))
    level,j=cells[i]
    assert level<12
    return sum(prove(cells[:i]+((level+1,v),)+cells[i+1:],reverse) for v in (2*j,2*j+1))


if __name__=="__main__":
    for reverse in (False,True):
        print("reciprocal coordinate",reverse,"certified boxes",prove(((0,0),)*3,reverse),
              "second / largest singular value above",target.str(4))
