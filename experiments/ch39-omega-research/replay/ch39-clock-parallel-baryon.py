"""A residual Z18 compatible with a clock-singlet Higgs and baryon triality on light singlets."""
from collections import defaultdict
from pathlib import Path
import runpy

z9=runpy.run_path(str(Path(__file__).with_name("ch39-clock-baryon-anomaly-free.py")))
m,sp=z9["m"],z9["sp"]

def check():
    fields=[("Q",6,2,3,1,1),("uc",3,1,0,-4,-1),("dc",3,1,0,2,-1),
            ("L",2,0,1,-3,0),("ec",1,0,0,6,0),("nc",1,0,0,0,0)]
    anomaly=[0]*6
    labels={}
    for name,dim,tc,tw,y,b in fields:
        h=(b+8*y)%9
        labels[name]=2*h
        src={0:{(0,0):3},3:m["A"],6:m["dual"](m["A"])}[h]
        p=m["shift"]({(2*a,2*w):v for (a,w),v in src.items()},y)
        spec=defaultdict(int)
        for (a,w),n in p.items():
            spec[a+4*w,w%18]+=n
        assert {k:v for k,v in spec.items() if v}=={(y,2*h):3}
        rank,a,_,a3,_,_,h3,ah2=m["moments"](p)
        assert rank==3
        anomaly=[v+u for v,u in zip(anomaly,[tc*a,tw*a,dim*a,dim*a3,dim*ah2,dim*h3])]
    assert anomaly[-2:]==[180144,-3888]
    correction=defaultdict(int)
    for target,mult in (((1296,1944),-139),((0,2916),94)):
        for key,v in z9["spectator"](6,target,mult,degree=9,step=72).items():
            correction[key]+=v
    spec=defaultdict(int)
    total=0
    for (p,q,a),v in correction.items():
        for h,n in sp["weights"](p,q).items():
            spec[a+4*h,h%18]+=v*n
            total+=abs(v)*n
            anomaly=[r+v*n*z for r,z in zip(anomaly,[0,0,a,a**3,a*h*h,h**3])]
    assert not any(spec.values()) and not any(anomaly)
    assert (3*labels["Q"]+labels["L"])%18==6
    print("Z18 light labels",labels,"UV anomaly totals",anomaly)
    print("additional spectator dimension with chirality multiplicity",total)
    return correction

if __name__=="__main__":
    check()
