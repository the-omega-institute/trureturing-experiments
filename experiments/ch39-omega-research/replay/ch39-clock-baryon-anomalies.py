"""Moments and continuous gauge anomalies of a concrete virtual SU(3)_X torus lift."""
from collections import defaultdict
from pathlib import Path
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-torus-baryon-lift.py")))
add=m["add"]
def conv(p,q):
    r=defaultdict(int)
    for (a,h),v in p.items():
        for (b,k),w in q.items():
            r[a+b,h+k]+=v*w
    return {k:v for k,v in r.items() if v}
one={(0,0):1}
fund={(0,1):2,(0,-2):1}
def dual(p):
    return {(-a,-h):v for (a,h),v in p.items()}
def shift(p,n):
    return {(a+n,h):v for (a,h),v in p.items()}
K=add((1,conv(fund,dual(fund))),(-5,one))
psi={(0,3*h):v for (_,h),v in fund.items()}
A=add((1,shift(conv(K,K),36)),(-8,shift(one,36)),(-2,shift(K,72)),
      (2,shift(psi,72)),(-1,shift(psi,36)))
A=shift(A,-12)
def moments(p):
    return [sum(v*g(a,h) for (a,h),v in p.items())
            for g in (lambda a,h:1,lambda a,h:a,lambda a,h:a**2,lambda a,h:a**3,
                      lambda a,h:h,lambda a,h:h**2,lambda a,h:h**3,lambda a,h:a*h**2)]
def check():
    # dim3, dim2, twice spectator color/weak indices (including multiplicity
    # of the other spectator), y=6Y, desired nontrivial discrete label.
    fields=[(6,2,3,1,True),(3,1,0,-4,False),(3,1,0,2,True),
            (2,0,1,-3,True),(1,0,0,6,True),(1,0,0,0,False)]
    anomalies=[0]*6
    for dim,tc,tw,y,has_charge in fields:
        p=shift(A if has_charge else {(0,0):3},y)
        rank,a,_,a3,h,h2,h3,ah2=moments(p)
        assert rank==3 and h==0
        vals=[tc*a,tw*a,dim*a,dim*a3,dim*ah2,dim*h3]
        anomalies=[n+v for n,v in zip(anomalies,vals)]
    print("single three-copy lift rank,A,A2,A3,H,H2,H3,AH2",moments(A))
    print("mixed color/weak A, gravitational A, A^3, AH^2, H^3",anomalies)
    return anomalies
if __name__=="__main__":
    check()
