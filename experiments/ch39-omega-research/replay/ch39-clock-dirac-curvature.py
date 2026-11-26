"""Minimum noncentral curvature and first heat correction forced by the two shell clocks."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
geo = runpy.run_path(str(root/"ch39-arc-relative-determinant.py"))
offset = arb(0,arb(".01"))
u = geo["a"]*geo["h0"]+offset
v = geo["b"]*geo["h0"]+offset
axis_sq = arb(8741547)/12563572
E = (1-axis_sq)*u.sin()**2*v.sin()**2
commutator = 1-2*E
even_distance,odd_distance = commutator.acos(),(-commutator).acos()
relative = 1+(odd_distance/arb.pi())**2
for name,z in (("even degree angle",even_distance),("odd degree angle",odd_distance),
               ("odd unit flux curvature heat coefficient ratio lower",relative.lower())):
    print(name,z.lower().str(12),z.upper().str(12))

# Uniform mass tower comparison after prescribing arbitrary SU(2) holonomies
# by a small determinant-one change of the Hermitian metric.
eps = arb(".0001")
assert arb(8).sqrt()*(-2*eps).exp()>arb(1).exp()
k = arb(2).sqrt()*(2*eps).exp()
s = arb("2.22")
assert 35*(-arb(".98")).exp()<s*(1+s**2)
assert k*(1+k**2*s**2)/(1+s**2)<arb(".96").exp()
print("max mass distortion after recalibrating the one fermion scale",(6*eps).exp().upper().str(12))

# The common-portal construction in the entire high-d region needs only
# a bound on total fermion weight. Recompute without an exact two-mode claim.
high=runpy.run_path(str(root/"ch39-common-portal-high-dirichlet.py"))
q=high["q"]
ferm_upper=q("36/5")*high["g2"].upper()**2*31
upper=(max(v.upper() for v in high["scalar_ratios"])
       +(ferm_upper+high["small"].upper())/high["den_lower"])/(8*high["D"])
b=high["lower"].lower().union(upper.upper())
f=high["f"]
kc,ks=b-f,b+3*f
margin=4*high["lh"]*q(".09")-kc**2
contact=2*109*q(".09")+(216*kc**2+72*ks**2).sqrt()
assert b>0 and b<f and margin>0 and contact<8*arb.pi()
print("curved fermion high-d common b",b.lower().str(12),b.upper().str(12),
      "contact",contact.upper().str(12))
