"""A weak-adjoint matter tower on the neutral graph and all-energy one-loop bounds."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q,r,sh = (short[k] for k in ("q","r","shell"))
pi = arb.pi()
k = 2*r*sum(short["total"])/pi
assert q("177.03844")<k<q("177.03845")
# Cancellable flat, Einstein and Ricci-flat a4 terms, per adjoint direction.
assert 3+1+4-2*(2+2) == 0
assert -F(1,2)+F(1,6)+4*(F(1,6)-F(1,4))+F(4,6) == 0
assert -F(1,15)+F(5,180)+F(28,720) == 0
assert -F(7,2)*2+F(2,6)+2*F(4,3) == -4
raw_extra = 3*(F(1,72)+F(1,72)+4*F(1,288))
assert raw_extra == F(1,8)

left, right = [], []
for _,_,d,_ in sh["edges"]:
    ell = r*(sh["rho"]*sh["theta"](d)).sqrt()
    base = pi/ell
    f = short["fsbox"][short["band"](d)]
    for n in range(1,6):
        left.append(n*base)
        right.append(n*base/f)
    assert 6*base/f>4 and 6*base>4
src = short["base"]["arc_budget"]()[1].upper()
neutral_top = (1+src/r**2).sqrt()
    
def upper_b(lo,hi,l,h):
    nl = 38+2*sum(not (m>hi) for m in left)
    dr = sum(not (m>hi) for m in right)
    nr = max(18,2*sum(m<lo for m in right))
    n0 = max(1,19 if neutral_top<lo else 1,
             sum(1+m**2<lo**2 for m in right))
    return q(41)/10+q(8)/3+l*nl+h*nr, q("18.5")+dr-4*n0

def check(l,h,ymax):
    # At x>=4 the elementary floor bounds are uniform and become stronger with x.
    assert l+h<0
    assert q(41)/10+q(8)/3+38*l-270*h+(l+h)*k*4<0
    weak_slope = k*(q(1)/2-q(15).sqrt()/2)
    assert weak_slope<0 and q("558.5")+4*weak_slope<0
    # Accumulate a bound on all possible losses of inverse alpha. This also
    # bounds every intermediate scale regardless of where a threshold lies.
    loss1=loss2=q(0)
    for j in range(192):
        lo=q(1)+q(j)/64
        hi=lo+q(1)/64
        b1,b2=upper_b(lo,hi,l,h)
        step=(hi/lo).log()/(2*pi)
        loss1+=max(q(0),b1.upper())*step
        loss2+=max(q(0),b2.upper())*step
    floor1=q("40.6142")-loss1
    floor2=q("40.6142")-loss2
    assert floor1>q(3)*ymax**2/5 and floor2>2
    print("uniform inverse alpha floors",floor1.lower().str(13),floor2.lower().str(13))
    print("asymptotic x alpha_1 and x alpha_2",(-2*pi/(k*(l+h))).str(12),(4*pi/(3*k)).str(12))
    return floor1, floor2

if __name__=="__main__":
    for lh,rh,ymax in ((144,4716,q(5)),(216,4068,q("4.5"))):
        check(q(lh)/90,q(rh-6048)/360,ymax)
    A0=k*pi.sqrt()/4
    print("additional pure R^2 UV slope",(A0*q(raw_extra)/(16*pi**2)).str(14))
