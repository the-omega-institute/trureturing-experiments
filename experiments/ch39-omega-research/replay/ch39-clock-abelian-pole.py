"""A necessary abelian running budget when the charge-two SU(3) completion turns on at Q."""
from pathlib import Path
import runpy
from flint import arb

short=runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q,r,sh=(short[k] for k in ("q","r","shell"))
pi=arb.pi()
gy2=12*pi/(5*q("40.9906"))
budget=8*pi**2*36/gy2
def levels(right):
    vals=[]
    for _,_,d,_ in sh["edges"]:
        length=r*(sh["rho"]*sh["theta"](d)).sqrt()
        if right:
            length*=short["fsbox"][short["band"](d)]
        fundamental=(pi/length).upper()
        assert fundamental>1
        for n in range(1,21):
            vals.append(n*fundamental)
        assert 21*fundamental>20
    return sorted(vals)

def exposure(x,upper_masses):
    # Each Dirichlet threshold is doubled; 18 modes are below Q on either side.
    return 18*x.log()+2*sum((x/v).log() for v in upper_masses[9:] if v<x)

def limit(generations,bR,base,left,right):
    lo,hi=q(1),q(20)
    def cost(x):
        return base*x.log()+generations*exposure(x,left)+bR*exposure(x,right)
    assert cost(hi)>budget
    for _ in range(40):
        mid=(lo+hi)/2
        if cost(mid)>budget:
            hi=mid
        else:
            lo=mid
    assert cost(hi)>budget
    return hi
def check():
    left,right=levels(False),levels(True)
    for label,bL,bR in (("four Y=1 Weyl",q(32)/3,q(2)),("six Y=1 Weyl",q(16),q(2)/3)):
        universal=limit(bL,bR,0,left,right)
        additional=limit(bL,bR,374,left,right)
        print(label,"vectorlike towers alone, then with minimal SM and two breaking fundamentals",
              universal.upper().str(12),additional.upper().str(12))
if __name__=="__main__":
    check()
