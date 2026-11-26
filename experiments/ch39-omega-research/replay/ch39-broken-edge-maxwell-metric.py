"""Intrinsic endpoint metrics after a symmetry-breaking positive single edge perturbation."""
from pathlib import Path
import runpy
from itertools import combinations
from flint import arb

work=Path(__file__).resolve().parent
edge=runpy.run_path(str(work/"ch39-clock-edge-source.py"))
rank=runpy.run_path(str(work/"ch39-clock-edge-rank-bound.py"))
shell=edge["s"]
src=edge["src"]
m=shell["matrix"](0,{})
J,W=src["J"],src["W"]
block=src["block"]
inv=block(m,J,J).inv()
base_sites=(11,12,18)
amp=block(m,base_sites,J)*inv
base=rank["c"][0,0][W.index(11),W.index(12)]
assert base.real>0
phase_error=arb(0)
for i in range(6):
    col=J.index(i)
    g=inv[col,col].real
    error=sum((abs(amp[j,col])*abs(amp[k,col])/(g*abs(base))).asin()
              for j,k in combinations(range(3),2))
    phase_error=max(phase_error,error.upper())
base_angle=3*base.arg()+arb(0,phase_error)
assert 0<base_angle<arb(".10822")

geo=runpy.run_path(str(work/"ch39-symmetric-shell-actions.py"))
density=[]
for f in src["faces"]:
    old=geo["geometry"]["lifted"](f,geo["negative"])
    phi=base_angle if set(f)==set(base_sites) else (arb.pi()-base_angle).union(arb.pi())
    density.append((phi-old)/geo["area"](f))
span=max(x.upper() for x in density)-min(x.lower() for x in density)
a=arb(18).sqrt()
b=a+2
s=b+7*(arb(2)/3).sqrt()
diameter=max(s.upper(),(arb(68).sqrt()+((s-a)**2+50).sqrt()).upper(),
             ((b**2+98).sqrt()+((s-b)**2+98).sqrt()).upper())
oscillation=diameter**2*span/4
cond=(arb(7930000)/4730000)/(1-arb(".148")**2)*(2*oscillation).exp()
ratio=arb(".095")/cond
print("base phase error, diameter and log-metric oscillation upper",
      phase_error.str(12),diameter.str(12),oscillation.upper().str(12))
print("new endpoint condition product upper and normalized Yukawa ratio lower",
      cond.upper().str(12),ratio.lower().str(12))
assert cond<arb("4.51") and ratio>arb(".021") and ratio>arb(".000271")/arb(".01392")
