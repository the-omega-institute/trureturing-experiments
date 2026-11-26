"""An integral tensor counterterm without charged spectators on the double cover of the determinant group."""
from pathlib import Path
from runpy import run_path
from fractions import Fraction as F

home=Path(__file__).resolve().parent
base=run_path(str(home/"ch39-shell-six-dimensional-anomaly.py"))
mon,basis,rank,compose=(base[k] for k in ("mon","basis","rank2","compose"))
# x = c1(det B) = 2y on the cover, y is primitive.
assert F(5,2)*4==10
# det B has degree 12, hence its root line has degree 6; the selected integral face lifts
# k_A=-3,-1 induce a -1 holonomy on this additional line for both kinds of faces.
turns=[3]*2+[5]*6
k=[-3]*2+[-1]*6
assert all((-v-j)%2==0 for v,j in zip(k,turns)) and -sum(k)//2==6
# One-polarization tensor with hyperbolic charge pairing: XY is integral.
# Purely four-dimensional gauge backgrounds extend trivially in the weak factors.
assert F(5,2)*24-3*9-4*8==1
assert 12//2==6
# On the cover Sq^2 c2=c3, Sq^2 c3=0 for A,B; y remains an independent degree-2 generator.
gens=[(mon(0,0),),(mon(5),),(mon(6),),(),(),(),(),()]
def cols(n):
    nxt=basis(n+1)
    out=[]
    for v in basis(n):
        terms=set()
        for j,t in enumerate(v):
            if t%2:
                for z in gens[j]:
                    terms.symmetric_difference_update([tuple(a+b-int(i==j) for i,(a,b) in enumerate(zip(v,z)))])
        out.append(sum(1<<nxt.index(t) for t in terms))
    return out
a,b=cols(2),cols(3)
assert rank(a)==2 and rank(b)==5 and len(b)==7 and not any(compose(i,b) for i in a)
print("index-two lift, extra determinant line degree 6, bare mixed anomaly integral Y coefficient",10)
print("fiber integral of Y = x4 = 2 y4; left weak anomaly in charge-one convention",6)
print("degree-seven spin bordism dimension",len(b)-rank(b)-rank(a))
# A faithful SU(39) trace metric: mixed, B times inverse root, inverse root, two adjoints.
branch=((12,1,3,4),(3,-2,0,1),(1,-6,0,0),(15,0,8,0),(8,0,0,6))
assert sum(d*q for d,q,_,_ in branch)==0
assert sum(d for d,_,_,_ in branch)==39
levels=(sum(a for _,_,a,_ in branch),sum(b for _,_,_,b in branch))
center=sum(F(d*q*q,144) for d,q,_,_ in branch)
assert levels==(11,11) and center==F(5,12) and F(1,36)<center/11<F(1,12)
print("faithful SU(39) trace levels and relative central stiffness",*levels,center/11)
