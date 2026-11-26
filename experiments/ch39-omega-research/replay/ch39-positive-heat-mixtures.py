"""Three species responses characterize a single exponential among positive heat mixtures."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root / "ch39-dense-portal-spectra.py"))
data = source["spectra"]()
q = source["q"]
y = data.mass_ratios
R2 = data.scale_lower.union(data.scale_upper)**2
gap = R2*(y[6]-y[0])
power = (y[18]-y[0])/(y[6]-y[0])
assert power>2
eps = q(".001")
# E under the measure tilted by exp(-t e_1); X = exp[-t(e_7-e_1)].
# Integrate a second-derivative Taylor remainder from the mean of X to 0.
# The coefficient is minimized at 0 since p>=2.
# Bound exp(2*p*eps)-1 divided by p-1 without losing the dependence on p.
# Use the anti-monotonicity of p/(p-1) to sharpen the dependency enclosure.
pl, ph = power.lower(), power.upper()
relative_variance = 2*eps*pl/(pl-1)*(2*eps*ph).exp()
for delta in (1,2,3):
    radius = 1-(-(gap.lower()*delta-2*eps)).exp()
    assert radius>0
    tail = relative_variance/radius**2
    print("tilted weight outside +/-",delta,"upper",tail.upper().str(12))
print("power",pl.str(12),ph.str(12),"relative variance upper",relative_variance.upper().str(12))
