"""A normalization-free species response predicted by the shortest-arc heat model."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root / "ch39-dense-portal-spectra.py"))
data = source["spectra"]()
q = source["q"]
R2 = data.scale_lower.union(data.scale_upper)**2
dy = data.mass_ratios[6]-data.mass_ratios[0]
gap = R2*dy
theta = lambda k: (1-q(k)/1734).acos()
time = theta(804)/theta(12)
relative_response = (time*gap).exp()
second_time = (16*gap).exp()
assert gap>0 and relative_response<second_time
# Perturbed ratios of the squared Dirichlet masses measured on finite
# segments versus continuum extremal edges. All h are scaled per edge.
h = q("1e-6")
lo = 1/(1+arb.pi()**2*h*h/6)
hi = 1/(1-h*h)**2
scale_change = (lo**2/hi).lower().union((hi**2/lo).upper())
assert abs(scale_change-1) < q("5.7e-12")
for name,v in (("normalized mass gap",gap),("arc-time response ratio",relative_response),
               ("time-16 response ratio",second_time),("upper log-error multiplier",1/gap)):
    print(name,v.lower().str(12),v.upper().str(12))
print("finite two-edge readout absolute deviation",abs(scale_change-1).upper().str(12))
