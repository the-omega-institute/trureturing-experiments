"""Trace-only rational-height selection including all sixteen cosine branches."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from flint import arb,acb,fmpz_mat
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-rational-height-selector.py")))
p,q,T,w,top,K=[m[k] for k in ("p","q","T","w","top","K")]
pi=arb.pi()
target=[(2*p*f).cos() for f in w]
read=[(2*q*T*f).cos() for f in w]
e0=sum((x-y)**2 for x,y in zip(read,target))
def inverse_width(v,eps):
    lo,hi=v-eps,v+eps
    assert -1<lo<hi<1
    a=v.acos()
    return max((a-hi.acos()).upper(),(lo.acos()-a).upper())/(2*pi)
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def round_to_int(v):
    return int((v+arb(1)/2).floor().unique_fmpz())
def certify(eps):
    widths=[inverse_width(v,eps) for v in target]
    assert K>w[0]*(q*top+p)/pi+widths[0]
    M=int((arb(K)/max(x.upper() for x in widths)).ceil().unique_fmpz())
    rr=[f/w[0] for f in w[1:]]
    A=[round_to_int(r*M) for r in rr]
    B=fmpz_mat([[1,*A],[0,M,0,0],[0,0,M,0],[0,0,0,M]])
    red,U=B.lll(transform=True)
    assert U*B==red and abs(U.det())==1
    basis=[[int(red[i,j]) for j in range(4)] for i in range(4)]
    stars=[]
    for v in basis:
        v=[F(z) for z in v]
        for u in stars:
            k=dot(v,u)/dot(u,u)
            v=[x-k*y for x,y in zip(v,u)]
        stars.append(v)
    short=min(dot(u,u) for u in stars)
    bound=(K**2+sum((M*(widths[0]*r+e)+arb(K)/2+1)**2 for r,e in zip(rr,widths[1:]))).upper()
    if not arb(short.numerator)/short.denominator>4*bound:
        return None
    for signs in product((-1,1),repeat=4):
        offset=[round_to_int(M*(signs[0]-s)*p*f/pi) for s,f in zip(signs[1:],w[1:])]
        target0=[0,*[-d for d in offset]]
        residual=target0[:]
        for v,u in zip(reversed(basis),reversed(stars)):
            k=round(dot(residual,u)/dot(u,u))
            residual=[x-k*y for x,y in zip(residual,v)]
        if arb(dot(residual,residual))<=bound:
            if signs!=(signs[0],)*4 or residual!=[0]*4:
                return None
    assert abs(T-arb(p)/q)<pi*widths[0]/(q*w[0])
    window=pi*widths[0]/(q*w[0])
    assert acb.zeta_zero(718).imag<arb(p)/q-window
    assert acb.zeta_zero(720).imag>arb(p)/q+window
    return bound
print("cosines at reference",*[v.str(10) for v in target])
print("719 total squared cosine deviation",e0.str(14))
for a in ("0.00001","0.000003","0.000001","0.0000003","0.0000001","0.00000003","0.00000001"):
    e=arb(a)
    if certify(e) is not None:
        print("all other positive zeros through 1e8 have max trace deviation above",a)
        assert e0<e*e
        break
else:
    raise ArithmeticError("no certified trace gap")
