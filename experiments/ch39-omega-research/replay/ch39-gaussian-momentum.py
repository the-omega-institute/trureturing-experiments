"""Continuum four-momentum average of the shell determinant at the eighth positive zero."""
from fractions import Fraction
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import runpy
from flint import arb, acb, acb_mat

root = Path(__file__).resolve().parent
shell = runpy.run_path(str(root/'ch39-arc-relative-determinant.py'))
def real(r):
    f=Fraction(r)
    return arb(f.numerator)/f.denominator
I=[6,7,8]
O=[i for i in range(19) if i not in I]
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
orbit=[6]
for _ in range(2):
    point=shell['sites'][orbit[-1]]
    orbit.append(shell['sites'].index(point[1:]+point[:1]))
fourier=[]
for s in range(3):
    vec=acb_mat(3,1)
    for j,i in enumerate(orbit):
        vec[I.index(i),0]=omega**(s*j)/arb(3).sqrt()
    fourier.append(vec)
g=acb.zeta_zero(8).imag
u=shell['a']*g+arb(0,'0.00001')
v=shell['b']*g+arb(0,'0.00001')
su,x=u.sin_cos()
sv,y=v.sin_cos()
@lru_cache(maxsize=256)
def edge_coeff(k,ds):
    if k==0:return real(1),real(1)
    t=(real(k)*shell['rho']*shell['theta'](ds)).sqrt()
    return t/t.sinh(),t/t.tanh()

def schur_terms(lo,hi):
    H=acb_mat(19,19)
    for i,V in enumerate(shell['V']):
        H[i,i]=V+real(lo).union(real(hi))
    weights={}
    W=[acb_mat(16,3) for _ in range(3)]
    for i,j,ds,z in shell['edges']:
        if ds not in weights:
            a,b=edge_coeff(lo,ds);c,d=edge_coeff(hi,ds)
            weights[ds]=a.union(c),b.union(d)
        off,diag=weights[ds]
        H[i,i]+=diag;H[j,j]+=diag
        H[i,j]=-off*z;H[j,i]=H[i,j].conjugate()
        if (i in I)!=(j in I):
            outer,inner=(i,j) if j in I else (j,i)
            assert ds not in (724,804) or i==outer
            W[{724:1,804:2}.get(ds,0)][O.index(outer),I.index(inner)]=H[outer,inner]
    A=acb_mat([[H[i,j] for j in O] for i in O])
    B=acb_mat([[H[i,j] for j in I] for i in I])
    inv=A.inv()
    sectors=[]
    for e in fourier:
        cols=[w*e for w in W]
        left=[v.conjugate().transpose() for v in cols]
        right=[inv*v for v in cols]
        const=(e.conjugate().transpose()*B*e)[0,0]
        const-=sum((l*r)[0,0] for l,r in zip(left,right))
        hs=[-(left[i]*right[j])[0,0] for i,j in ((0,1),(0,2),(1,2))]
        sectors.append((const.real,hs))
    return sectors

def log_derivative(sectors,c,order=1):
    z=x*y+c*su*sv
    result=real(0)
    for const,hs in sectors:
        a,b,d=(v.imag for v in hs)
        tr=const+2*sum(v.real*t for v,t in zip(hs,(x,y,z)))
        norm=(a*a*(1-x*x)+b*b*(1-y*y)+d*d*(1-z*z)
              +2*a*b*(z-x*y)+2*a*d*(x*z-y)+2*b*d*(x-y*z))
        det=tr*tr-4*norm
        if not det>0:raise ArithmeticError('sector determinant')
        deriv=4*tr*hs[2].real-8*(-d*d*z+a*b+a*d*x-b*d*y)
        first=su*sv*deriv/det
        if order==1:
            result+=first
        else:
            result+=(su*sv)**2*8*(hs[2].real**2+d**2)/det-first**2
    return result

def tail_bound(K):
    """Markov's inequality on the two-clock word and the unsigned adjacency walk counts."""
    edges={(i,j):ds for i,j,ds,_ in shell['edges']}
    mixed=0
    for i,j,k in combinations(range(19),3):
        ds=[edges.get(pair) for pair in ((i,j),(j,k),(i,k))]
        mixed+=int(None not in ds and 724 in ds and 804 in ds)
    adj=[{b if a==i else a for a,b in edges if i in (a,b)} for i in range(19)]
    trace4=sum(len(a&b)**2 for a in adj for b in adj)
    assert mixed==12 and trace4==47970 and max(map(len,adj))==18
    K=real(K)
    return 4*mixed/K+trace4/(K**2*(1-18/K))

def quadrature(mesh_error=None,order=1,tail_control=tail_bound):
    K=1024
    cs=[real(Fraction(i,8)).union(real(Fraction(i+1,8))) for i in range(-8,8)]
    total_lo=real(0)
    total_hi=real(0)
    core=real(0)
    loss=real(0)
    for lo,hi,count in [(0,1,64)]+[(2**i,2**(i+1),128 if i<4 else 384) for i in range(10)]:
        for j in range(count):
            a=Fraction(lo)+Fraction((hi-lo)*j,count)
            b=a+Fraction(hi-lo,count)
            terms=schur_terms(a,b)
            vals=[log_derivative(terms,c,order) for c in cs]
            if order==2:
                vals=[-v for v in vals]
            lower=min(v.lower() for v in vals)
            upper=max(v.upper() for v in vals)
            if mesh_error is not None:
                err=mesh_error(real(a),real(b)).upper()
                lower=(lower-err).lower()
                upper=(upper+err).upper()
            measure=real((b*b-a*a)/2)
            total_lo+=lower*measure
            total_hi+=upper*measure
            if hi<=16:
                assert lower>0
                core+=lower*measure
            else:
                loss+=min(real(0),lower)*measure
    tail=tail_control(K)
    for name,v in [('core to 16',core),('possible negative part through 1024',loss),
                   ('absolute integrated tail',tail),('all cutoffs beyond 16',core+loss-tail),
                   ('unregulated lower',total_lo-tail),('unregulated upper',total_hi+tail)]:
        print(name,v.lower().str(12),v.upper().str(12))
    assert core+loss>tail and total_lo>tail
    return (total_lo-tail).lower(),(total_hi+tail).upper()
if __name__=='__main__':
    quadrature()
