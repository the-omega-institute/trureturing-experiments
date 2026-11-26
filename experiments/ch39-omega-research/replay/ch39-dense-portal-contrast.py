"""Required variation of any positive common matrix, including noncommuting portals."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
bounds = runpy.run_path(str(root / "ch39-dense-common-portal-obstruction.py"))
q = bounds["q"]
U, eta, D = (bounds[k] for k in ("U","eta","D"))
base, linear, quadratic = (bounds[k] for k in ("base","linear","quadratic"))
n = bounds["nhi"]
allowed_ratio = q("2.15")
norm = 8*arb.pi()/q(48).sqrt()
for _ in range(6):
    average = (base+linear*norm+quadratic*norm**2)/(8*D)
    largest = allowed_ratio * average
    norm = (U/(6*eta)).sqrt()+q(n).sqrt()*largest
assert largest < q(".005")
print("contrast at most 2.15 forces maximum below",largest.upper().str(12))

# With heat time t, the inverse contrast is an average of exp(-t*(e_i-e_1)),
# which decreases at fixed spectrum. Use lower/upper spectral enclosures at a
# single t to cover the whole shorter-time interval.
data = bounds["spectrum"]["spectra"]()
ratio_sq = data.scale_lower.union(data.scale_upper)**2
ys = data.mass_ratios
energies = [y*ratio_sq for y in ys]
weights = [y*(1-y.log()) for y in ys]
time = q(".49")
response = sum(w.lower()*(-time*(v.upper()-energies[0].lower())).exp()
               for w,v in zip(weights,energies))
contrast = min(q(n),sum(weights).upper()) / response
assert contrast < allowed_ratio
print("heat time <= 0.49 has contrast below",contrast.upper().str(12))
