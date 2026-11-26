"""A nine-vertex stress certifies the 42-reading Euclidean obstruction."""
from fractions import Fraction as Q
from itertools import combinations,product
from math import prod
from pathlib import Path
from runpy import run_path
import json

from sympy import Matrix, Rational, zeros

directory=Path(__file__).parent
integer=run_path(str(directory/"ch39-rh-clock-integer-replay.py"))
cell,I=integer["read_cell"],integer["Bounds"]


def main():
    data=json.loads((directory/"ch39-rh-clock-integer-data.json").read_text())["inputs"]
    ref=list(map(cell,data["reference"]))
    pivot=cell(data["pivot"][0])
    weights=[cell(v)/(I(1)-pivot) for v in data["rho"]]
    nodes=[i*i for i in range(1,42)]
    products={s:prod(Q(s+v,s-v) for v in nodes if s!=v) for s in nodes}
    def response(t):
        exact=sum((w*c/(c+I(t)) for c,w in zip(ref,weights)),I(0))
        assert t not in nodes
        k=prod((t-s)/(t+s) for s in nodes)
        norm=sum(abs(2*s*k/(t-s)*products[s]) for s in nodes)
        err=2*abs(k)+norm/2**96
        return exact,I(exact.lo-err,exact.hi+err)
    pts=[p for p in product(range(5),repeat=3) if sum(x*x for x in p)==17]
    h,s,b=map(Q,("6.17986747834","21.97447633726","5.06908655386"))
    value={0:0,1:h,4:s}
    x=[[value[v]-(h+s)/3 for v in p] if len(set(p))==3 else
       [b*(int(v==3)-Q(1,3)) for v in p] for p in pts]
    edges={}
    for i,j in combinations(range(9),2):
        n=17+sum(a*b for a,b in zip(pts[i],pts[j]))
        kind=int(len(set(pts[i]))<3)+int(len(set(pts[j]))<3)
        edges.setdefault((kind,n),[]).append((i,j))
    active=((0,18),(0,33),(1,28),(1,31))
    matrices=[]
    pos=Matrix(x)
    for orbit in active:
        m=zeros(9)
        for i,j in edges[orbit]:
            m[i,i]+=1
            m[j,j]+=1
            m[i,j]-=1
            m[j,i]-=1
        matrices.append(m)
    coeff,=Matrix.hstack(*[(m*pos).reshape(27,1) for m in matrices]).nullspace()
    if coeff[0]>0:
        coeff=-coeff
    assert coeff[0]<0 and all(v>0 for v in list(coeff)[1:])
    stress=sum((v*m for v,m in zip(coeff,matrices)),zeros(9))
    assert stress*pos==zeros(9,3) and stress*Matrix([1]*9)==zeros(9,1)
    # det(lambda I + stress) has nonnegative coefficients. Symmetry makes
    # the eigenvalues real, so there can be no negative stress eigenvalue.
    assert all(v>=0 for v in (-stress).charpoly().all_coeffs()) and stress.rank()==6
    fixed_lo,fixed_hi=map(Q,(".06641861504",".06641861506"))
    all_lo,all_hi=map(Q,(".06641861",".06641863"))
    for orbit,pairs in edges.items():
        target=Q(orbit[1],34)
        for i,j in pairs:
            t=sum((a-b)**2 for a,b in zip(x[i],x[j]))
            exact,wide=response(t)
            for value,lo,hi in [(exact,fixed_lo,fixed_hi),(wide,all_lo,all_hi)]:
                assert target-hi<value.lo<=value.hi<target+hi
                if orbit==active[0]:
                    assert value.lo>target+lo
                elif orbit in active:
                    assert value.hi<target-lo
    print("Planar upper and all-dimensional lower error: (0.06641861504,0.06641861506).")
    print("Any positive response compatible with 96-bit data: (0.06641861,0.06641863).")


if __name__=="__main__":
    main()
