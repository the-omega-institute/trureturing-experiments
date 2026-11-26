"""Anomaly-free virtual fermion lift forbidding all single-lepton dimension-six proton operators."""
from pathlib import Path
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-baryon-anomalies.py")))
def check():
    fields=[("Q",6,2,3,1,0),("uc",3,1,0,-4,0),("dc",3,1,0,2,0),
            ("L",2,0,1,-3,1),("ec",1,0,0,6,-1),("nc",1,0,0,0,-1)]
    anomaly=[0]*6
    trial={}
    for name,dim,tc,tw,y,L in fields:
        p={0:{(0,0):3},1:m["A"],-1:m["dual"](m["A"])}[L]
        p=m["shift"](p,y)
        rank,a,_,a3,h,h2,h3,ah2=m["moments"](p)
        assert rank==3 and h==0
        restricted={}
        for (x,w),v in p.items():
            key=(x+4*w,w%9)
            restricted[key]=restricted.get(key,0)+v
        assert {k:v for k,v in restricted.items() if v}=={(y,(3*L)%9):3}
        trial[name]=(3*L)%9
        anomaly=[a+b for a,b in zip(anomaly,[tc*a,tw*a,dim*a,dim*a3,dim*ah2,dim*h3])]
    assert anomaly==[0]*6
    assert (3*trial["Q"]+trial["L"])%9 and (2*trial["uc"]+trial["dc"]+trial["ec"])%9
    assert (2*trial["uc"]+trial["dc"]+trial["ec"]+2*trial["nc"])%9==0
    assert 2*(-4)+2+6==0
    print("projected light charges",trial,"continuous perturbative anomalies",anomaly)
    print("dimension-six QQQL forbidden, dimension-nine uc uc dc ec nc nc allowed")
if __name__=="__main__":
    check()
