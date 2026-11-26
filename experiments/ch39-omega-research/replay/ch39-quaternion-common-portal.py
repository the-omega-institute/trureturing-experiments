"""One loop-order common portal in the eighteen-mode pole matching problem."""
import runpy
from pathlib import Path

from flint import arb

root = Path(__file__).resolve().parent
gauge = runpy.run_path(str(root / "ch39-quaternion-gauge.py"))
q = gauge["q"]


def interval(left, right):
    return q(left).union(q(right))


def endpoints(label, value):
    print(label, value.lower().str(15), value.upper().str(15))


# The 18-mode gauge certificate supplies this rectangle and the positive masses.
r = interval("31.0740", "31.6548")
scale = 1000 * interval("26.865", "26.873").exp()
alpha_inverse = interval("38.3183", "38.3273")
g2 = 4 * arb.pi() / alpha_inverse
exposures = gauge["geometry"](r)
mass_ratios = [x / r**2 for x in gauge["xs"]]
logs = [x.log() for x in mass_ratios]
assert len(logs) == 18 and all(t < -4 for t in logs)
W = -sum(logs) / 2
lh = q(".02335")
# The broad low-energy evolution contains the whole gauge rectangle.
ll = interval("-.031357", "-.0298109")
low_mass_sq = interval("10271.8", "10311.5")
D = 16 * arb.pi()**2
f = (2 * arb.pi()**2 * (lh - ll) / (3 * W)).sqrt()
fermion = q("36/5") * g2**2 * sum(
    (1 + 12*t) * (-2*t).exp() for t in exposures
)
scalar = sum(
    x * (q("78/25") * (3*l**2/2 + 2*l + q("7/2")) * g2**2
         + q("144/5") * (3-l) * g2*f + 24 * ((l-2)**2+1) * f**2)
    for x, l in zip(mass_ratios, logs)
)
hard = (scale/D)**2 * (scalar + fermion)
# b denotes the dimensionless common portal itself, counted at one loop.
# The derivative of the one-loop potential is linear in b at retained order.
response = 8 * scale**2 / D * sum(x*(1-l) for x, l in zip(mass_ratios, logs))
b = (hard + low_mass_sq) / response
kc, ks = b - f, b + 3*f
assert kc < 0 and ks > 0
color_quartic, neutral_quartic = q(".09"), q(".005")
stability = 4*lh*color_quartic - kc**2
# Three radial channels for 2 Higgs, 108 color and 36 neutral complex fields.
# The portal has off-diagonal radial norm sqrt(2*108*kc^2 + 2*36*ks^2).
radial = max(6*lh, 2*109*color_quartic, 2*37*neutral_quartic)
contact = radial + (216*kc**2 + 72*ks**2).sqrt()
assert max(abs(kc), abs(ks)) < 1 and stability > 0 and contact < 8*arb.pi()
assert hard > 0 and response > 0 and 0 < b and b < f
for label, value in (
    ("leading portal f", f), ("common portal b", b), ("unshifted mass GeV^2", hard),
    ("absolute slope GeV^2 per unit b", response),
    ("relative mass sensitivity", (hard+low_mass_sq)/low_mass_sq),
    ("stability margin", stability), ("tree contact upper", contact),
):
    endpoints(label, value)
