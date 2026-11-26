"""Exact matter cost of one Y=9/2 weak doublet in a layerwise semisimple completion."""
from fractions import Fraction as F
from math import prod

# The bound is measured in rank-one shell spectra of a pure YM adjoint per simple factor.
# A rank-two right matter layer appears twice as often asymptotically.
adj=(8,15,3,3,3)
ca=(3,4,2,2,2)
def xdim(p,q):
    return (p+1)*(q+1)*(p+q+2)//2
def cdim(u,v):
    return (u+1)*(v+1)*(u+2)*(v+2)*(u+v+3)//12
def examine(p,q,u,v,r,s):
    d=2*xdim(p,q)*cdim(u,v)*(r+1)*(s+1)
    weight=F(1,6) if p==q and u==v and (1+r+s)%2==0 else F(1,3)
    cas=(F(p*p+q*q+p*q+3*p+3*q,3),F(3*u*u+2*u*v+3*v*v+12*u+12*v,8),
         F(3,4),F(r*(r+2),4),F(s*(s+2),4))
    b=[weight*d*c/n for c,n in zip(cas,adj)]
    return max(6*x/(11*n) for x,n in zip(b,ca)),b
def target(p,q,u,v,r,s):
    hset={3*(i+j)-2*(p+2*q) for i in range(q,p+q+1) for j in range(q+1)}
    return any(abs(4*h+3*(v-u)+3*z+2*t)==27 for h in hset
               for z in range(-r,r+1,2) for t in range(-s,s+1,2))

def run():
    trial=(0,3,0,0,1,0)
    assert target(*trial)
    cut=examine(*trial)[0]
    # b_L >= (non-L dimension)/12, hence cost >= dimension/44.
    maxdim=int(44*cut)
    xs=[(p,q,xdim(p,q)) for p in range(maxdim+1) for q in range(maxdim+1)
        if xdim(p,q)<=maxdim]
    cs=[(u,v,cdim(u,v)) for u in range(maxdim+1) for v in range(maxdim+1)
        if cdim(u,v)<=maxdim]
    best=cut
    minima=[]
    for p,q,dx in xs:
        for u,v,dc in cs:
            lim=maxdim//(dx*dc)
            for r in range(lim):
                for s in range(lim//(r+1)):
                    rep=(p,q,u,v,r,s)
                    if not target(*rep):
                        continue
                    cost,b=examine(*rep)
                    if cost<best:
                        best,minima=cost,[]
                    if cost==best:
                        minima.append((rep,b))
    print("best common number of rank-one pure gauge towers",best,float(best))
    print("minimizing weights and per-matter-layer coefficients",minima)
    return best,minima
if __name__=="__main__":
    run()
