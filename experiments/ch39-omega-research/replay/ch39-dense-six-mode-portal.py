"""Solve the dense 172-mode problem with a common portal supported on the low spectral cluster."""
import runpy
from pathlib import Path

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root / "ch39-dense-portal-spectra.py"))
q = source["q"]
data = source["spectra"]()
assert data.count_lower == data.count_upper == 172
ys = data.mass_ratios
assert ys[5] < q(".046") and ys[6] > q(".092")
D = 16*arb.pi()**2
Q = 1000*q("26.845").union(q("26.889")).exp()
g2 = 4*arb.pi()/q(38).union(q("38.36"))
lh = q(".02335")
ll = q("-.031357").union(q("-.0298109"))
msq = q("10271.8").union(q("10311.5"))
ls = [y.log() for y in ys]
W = -sum(ls[:6])/2
f = (2*arb.pi()**2*(lh-ll)/(3*W)).sqrt()
pure = sum(y*q("78/25")*(3*l**2/2+2*l+q("7/2"))*g2**2 for y,l in zip(ys,ls))
leading = sum(y*(q("144/5")*g2*(3-l)*f+24*((l-2)**2+1)*f**2) for y,l in zip(ys[:6],ls[:6]))
# All allowed positive fermion exposures sum to one. Every t has w(t) >= 1,
# and w(t) <= 1+10t. At most 21 modes occur across the three Dirac lines.
fermion = q("36/5")*g2**2*q(1).union(q(31))
hard = (pure+leading+fermion)*(Q/D)**2
response = 8*Q**2/D*sum(y*(1-l) for y,l in zip(ys[:6],ls[:6]))
# Cancel the common Q^2 in the ratio before applying interval arithmetic.
b = (pure+leading+fermion+msq*(D/Q)**2)/(8*D*sum(y*(1-l) for y,l in zip(ys[:6],ls[:6])))
color, neutral = b-f, b+3*f
# Each pure heavy quartic couples the radial variables only within its group.
# The 166 remaining levels have zero Higgs portal.
color_self = q(".18")
neutral_self = q(".005")
rest_self = q(".001")
margin = 4*lh*color_self - max(q(0), (-color).upper())**2
bare = max(6*lh,2*37*color_self,2*13*neutral_self,2*(166*8+1)*rest_self)
contact = bare+(72*color**2+24*neutral**2).sqrt()
for name,z in (("f",f),("b",b),("color",color),("neutral",neutral),
               ("stability margin",margin),("contact",contact),
               ("mass sensitivity",1+hard/msq),("unexpanded quartic relative term",b**2/(3*f**2))):
    print(name,z.lower().str(12),z.upper().str(12))
assert b>0 and neutral>0 and margin>0 and contact<8*arb.pi()

# A second fixed shape couples the other 166 levels at one thirty-second
# of the strength on the six lowest modes. The leading F is still f P_6.
eps = q("1/32")
low_response = sum(y*(1-l) for y,l in zip(ys[:6],ls[:6]))
rest_response = sum(y*(1-l) for y,l in zip(ys[6:],ls[6:]))
b = (pure+leading+fermion+msq*(D/Q)**2)/(8*D*(low_response+eps*rest_response))
color, neutral, rest = b-f, b+3*f, eps*b
color_self = q(".333")
margin = 4*lh*color_self-color**2
radial_diags = [74*color_self, 26*neutral_self, 2*(166*8+1)*rest_self]
coupling_squares = [72*color**2, 24*neutral**2, 2*166*8*rest**2]
cap = q("24.76")
for sign in (-1,1):
    assert all(cap+sign*x > 0 for x in radial_diags)
    assert cap+sign*6*lh > sum(w/(cap+sign*x) for w,x in zip(coupling_squares,radial_diags))
# The traceless and mixed-species channels have magnitude bounded by the
# individual portal couplings and the intra-species non-radial contact 2*c.
assert max(abs(color),abs(neutral),abs(rest),2*color_self) < 1
assert cap < 8*arb.pi() and color<0 and rest>0 and neutral>0 and margin>0
relative_quartic = (b/f)**2*(1+eps**2*sum(ls[6:])/sum(ls[:6]))/3
for name,z in (("1/32 shape common b",b),("actual color stability margin",margin),
               ("unexpanded quartic relative term",relative_quartic)):
    print(name,z.lower().str(12),z.upper().str(12))
