"""Two nonzero neutral thresholds with the mass ratio fixed by shell alignment."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

base=runpy.run_path(str(Path(__file__).with_name("ch39-three-threshold-newton.py")))
q,pi,A,Q,mp,baseline=[base[x] for x in ("q","pi","A","scale","mp","baseline")]
couplings=(F(1,6),F(2,5),F(13,60))
assert [F(1,2)-2*x for x in couplings]==[F(1,6),-F(3,10),F(1,15)]
L=48*pi.sqrt()
B=36*(q(9)/4).log()
def moment(nu):
    return A*L*nu*nu.sqrt()+(-q(135)/2).union(q(19))*B*nu
target=60*pi*((mp/Q)**2-baseline)
def bisect(is_low):
    lo,hi=q(0),q(10)**6
    assert moment(hi).lower()>target.upper()
    for _ in range(60):
        mid=(lo+hi)/2
        passed=moment(mid).upper()<target.lower() if is_low else moment(mid).lower()<=target.upper()
        if passed:
            lo=mid
        else:
            hi=mid
    return lo if is_low else hi
lo,hi=bisect(True),bisect(False)
assert moment(lo).upper()<target.lower() and moment(hi).lower()>target.upper()
nu=lo.lower().union(hi.upper())
unit_nu=198/q(17*867).sqrt()
relative=(Q/mp)**2*(baseline+moment(unit_nu)/(60*pi))
ratio=nu/unit_nu
for l,z in (("matching nu",nu),("first mass/Q",2*nu.sqrt()),("last mass/Q",3*nu.sqrt()),
            ("dimensionless stiffness over inertia needed",ratio),
            ("fraction of inverse Newton constant at unit stiffness and inertia",relative)):
    print(l,z.lower().str(13),z.upper().str(13))
