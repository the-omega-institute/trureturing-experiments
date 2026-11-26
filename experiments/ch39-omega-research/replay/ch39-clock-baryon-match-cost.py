"""Perturbative matching cost of the baryon-protecting clock-singlet Higgs and chiral lifts."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-baryon-anomalies.py")))
pi=arb.pi()
gY2=12*pi/(5*arb("40.9906"))
# For Y=(4H+a)/6, squared norm of H in SU(3) generator normalization is 12.
assert F(4,9)*12==F(16,3)
higgs_h,higgs_a,dim=18,-72,55
higgs_index=F(495,2)
bH=F(dim*higgs_a**2,3)
r,s2=m["moments"](m["A"])[0],4*m["moments"](m["A"])[2]
assert r==3 and s2==-17280
min_a2=3*abs(2*s2+r*(16+4))+2*abs(s2+r*9)
bmin=F(2*min_a2,3)+bH
assert min_a2==138006 and bmin==187044
# Abelian charge is normalized as the integer a in the covariant derivative.
charge=higgs_a**2*gY2/36
loop=arb(int(bmin))*gY2/(36*16*pi**2)
span=8*pi**2*36/(arb(int(bmin))*gY2)
# Isolated scalar-five-and-fifty contribution alone is unavoidable at a common threshold.
universal=8*pi**2*36/(arb(int(bH))*gY2)
print("square of physical Higgs abelian coupling lower bound",charge.lower().str(13))
print("loop parameter lower, one-loop pole E / Q upper", loop.lower().str(13),span.exp().upper().str(13))
print("Higgs alone: bA and common-threshold pole E/Q upper",bH,universal.exp().upper().str(13))
