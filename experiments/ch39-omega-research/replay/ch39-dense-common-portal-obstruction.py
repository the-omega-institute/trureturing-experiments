"""Incompatible mass and quartic bounds for an arbitrary leading F and a uniform common b."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
spectrum = runpy.run_path(str(root / "ch39-dense-portal-spectra.py"))
q = spectrum["q"]
rlo, rhi, nlo, nhi, log_bounds, common_response, _ = spectrum["spectra"]()
assert nlo == nhi == 172 and common_response > 150
D = 16*arb.pi()**2
g2 = 4*arb.pi()/38
lh = q(".02335")
U = 16*arb.pi()*lh
eta = 1 + 4*lh/(8*arb.pi()-6*lh)
target = 2*arb.pi()**2*(lh+q(".0298109"))/3
Qmin = 1000*q("26.845").exp()

# Upper coefficients in D^2 Delta_0 / (Q^2 sum y_i*(1-log y_i)).
pure = []
for j in range(512):
    l = (-q("3.58")*j/512).union(-q("3.58")*(j+1)/512)
    pure.append((q("78/25")*(3*l**2/2+2*l+q("7/2"))*g2**2/(1-l)).upper())
assert max(log_bounds) < q("3.58")/2
base = max(pure) + (q("36/5")*g2**2*31 + q("10311.5")*(D/Qmin)**2)/common_response
linear = 3*q("144/5")*g2/common_response.sqrt()
quadratic = 24*(5+2/q(1).exp()+arb.pi()**2/6)/common_response

# The squared row norm from the Higgs radial state to the heavy states
# is 2*(24*||F||_HS^2 + 8*n*b^2). Each iteration gives a necessary bound.
norm = 8*arb.pi()/q(48).sqrt()
for _ in range(6):
    b_upper = (base+linear*norm+quadratic*norm**2)/(8*D)
    norm = (U/(6*eta)).sqrt()+q(nhi).sqrt()*b_upper
assert b_upper < q(".003")

# The negative eigenvalues of the actual portals imply
# sum (6*eta*v^2+6*v^4/U) < U, v=(|f|-b_+)_+.
# In the mass basis, concavity of log bounds T from above by Tr(L F^2).
# Rearrangement permits matching the largest L and |f|. A dual upper bound
# only needs a scalar quartic maximum for each resulting L.
def quartic_bound(b, multiplier):
    total = multiplier*U
    k = 6*eta*multiplier
    c = 6*multiplier/U
    for enclosed in log_bounds:
        L = enclosed.upper()
        def derivative(v):
            return 2*(L*b+(L-k)*v-2*c*v**3)
        lo, hi = q(0), q(2)
        assert derivative(lo)>0 and derivative(hi)<0
        # The derivative first increases then decreases, or decreases throughout.
        # It has exactly one positive zero, which is the global maximizer.
        for _ in range(46):
            mid = (lo+hi)/2
            sign = derivative(mid)
            if sign>0:
                lo = mid
            else:
                assert sign<0
                hi = mid
        v = lo.union(hi)
        total += L*(v+b)**2-k*v**2-c*v**4
    return total

ceiling = quartic_bound(q(".005"), q(".279"))
assert ceiling < target
print("gauge scale ratio", rlo.lower().str(10), rhi.upper().str(10), "active", nlo, nhi)
print("mass requires b below", b_upper.upper().str(14))
print("b<=0.005 quartic maximum", ceiling.upper().str(14), "target floor",target.lower().str(14))
