"""A point-inertia family gauge bound independent of the relative masses of the two shells."""
from pathlib import Path
from runpy import run_path
from fractions import Fraction as F
from flint import arb

m=run_path(str(Path(__file__).with_name("ch39-shell-point-inertia-gravity.py")))
q, A, Q, mp, baseline=[m[k] for k in ("q","A","Q","mp","baseline")]
pi=arb.pi()
ia,it=m["moments"](17)
ja,jt=m["moments"](867)
rho=q(F(1601,1032))
assert rho>1
# All squared mass shifts scaled by Q^2. For fixed larger mass a,
# b <= a/(28/27), and the exact moment increases in both arguments.
u=q(F(27,28))
J=A*(4*pi.sqrt()/3)*((1+u)**(q(3)/2)-1-u**(q(3)/2))
L=19*((1+u)*(1+u).log()-u*u.log())


def upper(a):
    return (Q/mp)**2*(baseline+(J*a*a.sqrt()+L*a)/(12*pi))


lo,hi=q(0),q(10)**7
for _ in range(65):
    mid=(lo+hi)/2
    if upper(mid).upper()<1:
        lo=mid
    else:
        hi=mid
assert upper(lo).upper()<1
cap=16*pi/(q(it/ia).sqrt()+q(jt/ja).sqrt())**2
gmax=(cap/lo).sqrt()
gu=(4*pi/q("40.9906")).sqrt()
print("larger squared offset lower bound",lo.str(14),"lighter mass /Q necessary",(lo/rho).sqrt().str(13))
print("axial squared offset upper times g_F^2 if all gauge masses <=Q",cap.str(14))
print("necessary g_F upper",gmax.str(12),"two family vectors at g_U heavier than Q times",(gu/gmax).str(12))
for ratio in (q(1),q(9)/495):
    g2=gu**2*ratio
    frac=upper(cap/g2)
    print("family coupling squared / unified coupling squared",ratio.str(8),
          "inverse Newton upper fraction for vectors below Q",frac.upper().str(12))
    assert frac.upper()<1
