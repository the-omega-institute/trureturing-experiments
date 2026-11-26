"""Arithmetic distribution and source strength for prescribed rank-two torus clocks."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
geo = runpy.run_path(str(root / "ch39-arc-relative-determinant.py"))
pi = arb.pi()
c2 = arb(8741547)/12563572
c = c2.sqrt()
s = (1-c2).sqrt()
dmin = 2*c.asin()
floor = 1+(dmin/pi)**2
phase_pad = arb(0,arb(".01"))
p = geo["a"]*geo["h0"]+phase_pad
q = geo["b"]*geo["h0"]+phase_pad
v = p.sin()**2
w = q.sin()**2
# Stable expression for the smaller eigenvalue of the two-axis Gram bound.
mu = 2*(1-c2)*v*w/(v+w+((v-w)**2+4*c2*v*w).sqrt())
dodd = (2*(1-c2)*v*w-1).acos()
source_floor = dodd/(1/(2*mu).sqrt()+1/arb(2).sqrt())
print("global infimal action / central action",floor.lower().str(13),floor.upper().str(13))
print("area^(3/2) source lower in zero-35 square",source_floor.lower().str(13))
print("lower-end density coefficient",(pi*c/(2*dmin*s)).str(13),
      "upper-end log coefficient",(1/(pi*s)).str(13))

def probability(ratio, steps=4096):
    cutoff = (pi*(ratio-1).sqrt()/2).cos()/s
    assert 0<cutoff and cutoff<1
    start = cutoff.asin()
    dx = (pi/2-start)/steps
    # Monotone left/right sums. Endpoint values are known exactly.
    interior = sum((cutoff/(start+i*dx).sin()).acos() for i in range(1,steps))
    low = 4*dx*interior/pi**2
    high = low + 4*dx*cutoff.acos()/pi**2
    return low.lower().union(high.upper())

for z in ("1.5","1.55","1.6"):
    pr = probability(arb(z))
    print("positive-zero fraction, ratio <=",z,pr.lower().str(12),pr.upper().str(12))
