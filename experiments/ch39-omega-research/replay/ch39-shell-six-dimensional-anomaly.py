"""Six-dimensional bulk anomaly cancellation and the spin-bordism check for the connected product group."""
from fractions import Fraction as F
from itertools import product

# Local anomaly [Ahat (ch(W)ch(L)-ch(W*)ch(R))]_8. All weak-independent
# terms cancel. Rank and ch2 of W determine the remainder.
# Coefficient order: rank, ch2(A), ch2(B), x^2, where x=c1(B)=-c1(A).
blocks={
    "mixed":(12,3,4,-1),
    "A":(4,1,0,0),
    "B":(3,0,1,0),
    "one":(1,0,0,0),
    "detB":(1,0,0,F(1,2)),
    "dual_detB":(1,0,0,F(1,2))}
coeff={"mixed":1,"A":-3,"B":-4,"one":10,"detB":1,"dual_detB":1}
assert [sum(coeff[key]*block[j] for key,block in blocks.items()) for j in range(4)]==[0]*4
charge=dict(mixed=1,A=-3,B=4,one=0,detB=12,dual_detB=-12)
copies={r:coeff[r]*block[0]*charge[r] for r,block in blocks.items()}
assert sum(abs(v) for r,v in copies.items() if r!="mixed")==108
assert sum(v*charge[r] for r,v in copies.items())==0
print("local anomaly rank and second Chern character zero; extra absolute index per weak side",108)

# Compute the dual of the spin AHSS d2 at the sole total-degree-seven candidate.
# Integral cohomology is polynomial; reduction in degree eight is surjective.
# Chern generators x,a,b,l,r,A3,B3,A4 have degrees in complex units below.
weights=(1,2,2,2,2,3,3,4)
def basis(n):
    return [v for v in product(*(range(n//k+1) for k in weights)) if sum(a*b for a,b in zip(v,weights))==n]
def mon(*indices):
    return tuple(indices.count(i) for i in range(len(weights)))
# Sq^2 c_j=c1*c_j+(j-1)c_{j+1} mod 2 for the degrees in use.
gens=[(mon(0,0),),(mon(0,1),mon(5)),(mon(0,2),mon(6)),(),(),(mon(0,5),),(mon(0,6),),()]
def square(v):
    out=set()
    for j,mul in enumerate(v):
        if mul%2:
            for z in gens[j]:
                term=tuple(x+y-int(i==j) for i,(x,y) in enumerate(zip(v,z)))
                out.symmetric_difference_update([term])
    return out
def columns(n):
    b=basis(n+1)
    return [sum(1<<b.index(p) for p in square(v)) for v in basis(n)]
def rank2(cols):
    pivots={}
    for x in cols:
        while x:
            k=x.bit_length()-1
            if k not in pivots:
                pivots[k]=x
                break
            x^=pivots[k]
    return len(pivots)
def compose(vec,cols):
    out=0
    for i,z in enumerate(cols):
        if vec & (1<<i):
            out^=z
    return out
a,b=columns(2),columns(3)
assert all(compose(x,b)==0 for x in a)
assert (len(a),len(b),rank2(a),rank2(b))==(5,7,2,5)
print("cohomology ranks for Sq^2 from degree 4 and 6",rank2(a),rank2(b),
      "remaining spin bordism degree seven",len(b)-rank2(a)-rank2(b))
