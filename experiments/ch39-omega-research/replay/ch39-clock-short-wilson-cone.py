"""A simultaneous separation certificate against all positive adjoint loops of length at most four."""
from functools import lru_cache
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

loop=runpy.run_path(str(Path(__file__).with_name("ch39-clock-wilson-length.py")))
m,q=loop["m"],loop["q"]
pi=arb.pi()
half=q(".0005")
angles=[(q(0),q(0)),(pi/2,37*pi/24)]

def rotations(p,t):
    (s,c),(v,b) = p.sin_cos(),t.sin_cos()
    h,k = m["m"]["c"],m["m"]["v"]
    return {724:acb_mat([[acb(c,s),0],[0,acb(c,-s)]]),
            804:acb_mat([[acb(b,h*v),k*v],[-k*v,acb(b,-h*v)]])}

@lru_cache(None)
def prices(word):
    vals=[]
    pad=arb(0,half.upper())
    ref=m["rotations"](719)[0]
    for T in [ref]+[rotations(p+pad,t+pad) for p,t in angles]:
        U=loop["I"]
        for i in word:
            mat=T[(724,804)[abs(i)-1]]
            U=U*(mat if i>0 else mat.conjugate().transpose())
        vals.append((U[0,0]+U[1,1]).real)
    return vals

def check():
    gains=[]
    for n in (3,4):
        for _,_,z,word in loop["cycles"](n):
            if word:
                a,b,c=[8-t*t-2*z*t for t in prices(word)]
                gains.append((a-(b+2*c)/3).lower())
    print("all active 3/4-cycle margin",min(gains).str(12),"cycles",len(gains))
    assert min(gains)>q(".25")
    print("individual square zero density",(half*half/pi**2).lower().str(12))
if __name__=="__main__":
    check()
