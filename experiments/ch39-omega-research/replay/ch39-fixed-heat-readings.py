"""Fixed heat-portal parameters at the 35th zero. Bound tolerances at recurrent heights."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
heat = runpy.run_path(str(root / "ch39-dense-heat-portal.py"))
q = heat["q"]
data = heat["data"]
d = q(".5")
pi = arb.pi()
theta = lambda k: (1-q(k)/1734).acos()
t = theta(804)/theta(12)
f,b,*_ = heat["heat_match"](t)
R2 = data.scale_lower.union(data.scale_upper)**2
Q = heat["Q"]
g2 = heat["g2"]
D = heat["D"]
h = q("1e-6")
# The covariant difference per unit phase is bounded by the trapezoid
# mass on clock edges. Their maximum degree is four, with endpoint half-cells.
arc_norm = 1/(pi**2*theta(724)/theta(804)) + 2*h/d
assert arc_norm < q(".12") and data.scale_upper**2 < q("2.1")
phase_rate = 2*(arc_norm*(q("2.1")+10/d)).sqrt()+arc_norm*q(".01")
y_rate = phase_rate / data.scale_lower**2

quartic_rate = arb(0)
hard_rate = arb(0)
W_rate = arb(0)
for box in data.mass_ratios:
    # A one-dimensional interval contains the interpolating y between
    # either value at fixed Q,M,R,g,f,b. No moving scalar threshold occurs.
    l = box.log()
    p = (-t*R2*box).exp()
    dp = -t*R2*p
    portal = f*p
    dportal = f*dp
    A = q("78/25")*(3*l**2/2+2*l+q("7/2"))
    B = q("144/5")*(3-l)
    C = 24*((l-2)**2+1)
    # d[y A(log y)]/dy = A + dA/dl, likewise for the other polynomial terms.
    deriv = ((A+q("78/25")*(3*l+2))*g2**2
             + (B-q("144/5"))*g2*portal + box*B*g2*dportal
             + (C+48*(l-2))*portal**2 + 2*box*C*portal*dportal
             + 8*D*b*(p*l+box*dp*(l-1)))
    hard_rate += abs(deriv).upper()
    quartic_rate += abs(f**2*(p**2/box+2*l*p*dp)).upper()
    W_rate += (1/(2*box)).upper()
hard_rate *= (Q/D)**2*y_rate
quartic_rate *= 3*y_rate/(4*pi**2)
W_rate *= y_rate
delta = q("1e-29")
print("mass change GeV^2 per radian",hard_rate.upper().str(14))
print("absolute quartic per radian",quartic_rate.upper().str(14))
print("scalar exposure per radian",W_rate.upper().str(14))
while hard_rate*delta >= 1:
    delta /= 10
assert hard_rate*delta<1 and quartic_rate*delta<q("1e-15") and W_rate*delta<q("1e-12")
print("admissible phase tolerance",delta.str(10),"lower zero density",(2*(delta/pi)**2).lower().str(12))
