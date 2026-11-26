"""Eighteen active levels for every clock pair at high first Dirichlet energy."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / "ch39-quaternion-gauge.py"))
q = model["q"]
pi = arb.pi()
D = 16*pi**2


def interval(left, right):
    return q(left).union(q(right))


# Coarse endpoints from both all-clock gauge enclosures.
W = interval(47, "49.2")
g2 = 4*pi / interval(38, "38.36")
lh = q(".02335")
ll = interval("-.031357", "-.0298109")
msq = interval("10271.8", "10311.5")
Q = 1000*interval("26.845", "26.889").exp()
d0 = q(1225)
energy = q("34.3")**2
k = pi*(energy/d0).sqrt()
shell = runpy.run_path(str(root / "ch39-arc-relative-determinant.py"))
assert len(shell["V"]) == 19 and min(shell["V"][9:]) >= 959990
assert 959990 - energy - 36*k/k.sin() > 0
beta = pi**2 / (6*100**2*d0)
assert energy/(1+beta*(energy+10)) > q("34.181")**2
# On the low 18-dimensional trial space each block diagonal is at most 8
# and each site has at most 18 neighbors with unit coupling in operator norm.
assert q(".919").sqrt() * (q(47)/18).exp() > 13
assert (q(13)/q(26).sqrt()).log() > q(".5")

# The two shortest triangles have exactly two active levels throughout this
# r range. The geometry routine also tests the next level of every line.
weights = []
for j in range(88):
    r = (q(13)+q(j)/4).union(q(13)+q(j+1)/4)
    t = model["geometry"](r)
    weights.append(sum((1+12*v)*(-2*v).exp() for v in t))
weight = min(v.lower() for v in weights).union(max(v.upper() for v in weights))

f = (2*pi**2*(lh-ll)/(3*W)).sqrt()
# The mean of y*(1-log(y)), y = exp(-2*t), has a positive Jensen lower
# bound since each scalar t > 1/2 and mean(t) = W/18.
mean_t = q("49.2")/18
den_lower = 18 * (-2*mean_t).exp() * (1+2*mean_t)
assert den_lower > q(".4918")
scalar_ratios = []
for j in range(1024):
    l = (-q("7.18")*j/1024).union(-q("7.18")*(j+1)/1024)
    A = q("78/25")*(3*l**2/2+2*l+q("7/2"))
    B = q("144/5")*(3-l)
    C = 24*((l-2)**2+1)
    scalar_ratios.append((A*g2**2+B*g2*f+C*f**2)/(1-l))
ferm = q("36/5")*g2**2*weight
small = msq*(D/Q)**2
upper = (max(v.upper() for v in scalar_ratios) + (ferm+small).upper()/den_lower)/(8*D)
# Each scalar term is strictly positive: 1.5*l^2+2*l+3.5 >= 17/6.
# The fermion term alone bounds the affine cancellation sensitivity.
lower = min(v.lower() for v in scalar_ratios)/(8*D)
b = lower.lower().union(upper.upper())
assert b > 0 and b < f
kc, ks = b-f, b+3*f
margin = 4*lh*q(".09") - kc**2
contact = 2*109*q(".09") + (216*kc**2+72*ks**2).sqrt()
relative = 1 + ferm*(Q/D)**2/msq
assert margin > 0 and contact < 8*pi
for name, val in (("f", f), ("b", b), ("stability", margin),
                  ("contact upper", contact.upper()), ("sensitivity lower", relative.lower())):
    print(name, val.lower().str(12), val.upper().str(12))
