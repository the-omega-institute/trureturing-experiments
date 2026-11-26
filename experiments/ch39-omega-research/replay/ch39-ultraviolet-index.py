"""Necessary one-loop adjoint index for an unbroken simple group above the r=400 crossing."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-slope-matched-one-loop.py"))
q=m["q"]
pi=arb.pi()
h=m["h"]

def scalar_logs(ratio, count):
    summation=arb(0)
    for mult,edge_d in m["factors"]:
        for k in range(1,count+1):
            mr=k*edge_d.sqrt()/(400*(1-k*k*h*h))
            assert mr>1
            if mr<ratio:
                summation+=mult*(ratio/mr).log()
            else:
                assert mr>ratio
    return summation

for scale in (5,6,10,10000):
    p=q(scale)
    total=scalar_logs(p, min(99,scale))
    # Necessary even just for positive inverse coupling, irrespective of index I.
    minimum=(2*total-6*pi*q("33.2625"))/(11*p.log())
    print("P/Q",scale,"sum logs at least",total.lower().str(12),"I*C_A above",minimum.lower().str(12))
    assert minimum > q({5:"18.7537", 6:"32.4018", 10:"75.4816", 10000:"2391.3596"}[scale])
