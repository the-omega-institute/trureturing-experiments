"""Perturbative horizon above the trace-matched r=400 crossing if the product theory is continued."""
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-slope-matched-one-loop.py"))
q=m["q"]
pi=arb.pi()
r=q(400)
# Boundary-zero modes on distinct edges combine in one trial subspace.
# Below the candidate scale at most four harmonics per edge need be used.
def weak_upper(ratio):
    budget=q(0)
    for multiplicity,edge_d in m["factors"]:
        for j in range(1,5):
            relative=j*edge_d.sqrt()/(r*(1-j*j*m["h"]**2))
            assert relative>1
            if relative<ratio:
                budget+=multiplicity*(ratio/relative).log()
            else:
                assert relative>ratio
    yukawa_bound=3*4*pi*q("7/2")/(16*pi**2)
    alpha=q("33.2625")+((q("19/6")+yukawa_bound)*ratio.log()
                          -q("2/3")*budget)/(2*pi)
    return budget,alpha
for endpoint in ("3.9","4"):
    E,alpha=weak_upper(q(endpoint))
    print(endpoint,"additional scalar exposure lower",E.lower().str(13),"weak inverse upper",alpha.upper().str(13))
assert weak_upper(q("3.9"))[1]<1
