"""Finite weights for the four pure clock orbits and CP cancellation by sector exchange."""
from pathlib import Path
import runpy
from flint import arb

phase=runpy.run_path(str(Path(__file__).with_name("ch39-symmetric-source-transition.py")))
B=phase["B"]
mats=[phase["data"](i) for i in range(4)]
best,_,coeff=mats[3]
ratios=list(map(arb,("1.0004629","1.0001020","1.0008197")))
for (det,_,vals),ratio in zip(mats[:3],ratios):
    a,c,d=[det/best*x-ratio*y for x,y in zip(vals,coeff)]
    assert a+min(arb(0),c.lower())*B**2/3+min(arb(0),d.lower())*(B/3)**3>0
p=1/(1+sum(r**(-30000) for r in ratios))
assert p>arb(".955")
print("lower probability of the 804 orbit among four equiprobable pure choices, 30000 copies",p.str(12))
