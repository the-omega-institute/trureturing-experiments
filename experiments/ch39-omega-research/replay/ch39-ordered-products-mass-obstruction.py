"""Arbitrary products of small one-step defects about the tenth-zero ordered clock background."""
from pathlib import Path
import runpy
from flint import arb

phase=runpy.run_path(str(Path(__file__).with_name("ch39-symmetric-source-transition.py")))
fit=phase["fit"]
_,k,_=phase["data"](3)
t=-fit["block"](k,fit["ii"],fit["oo"])
h,e1,e2=[abs(v) for v in fit["symbols"](t)]
assert h>e1>e2
dist=arb("1.809")
Y=-k[fit["W"].index(11),fit["W"].index(12)]
assert arb(".034")<Y.arg()<arb(".038")
assert phase["finite"]["flux"]["degrees"](k)[0]==(3,[3])


def test(radius,slope,n,middle_floor,light_ceiling):
    eta=arb(radius)
    r=arb(slope)
    norm=(1+r*r).sqrt()
    growth=h.lower()-eta*norm
    assert e1*r+eta*norm<r*growth
    # A two-dimensional volume grows by at most this factor per step.
    upper=(h.upper()+eta)*(e1.upper()+eta)/growth**2
    lower=(e2.lower()-eta)*growth/((h.upper()+eta)*(e1.upper()+eta))
    assert 0<lower<1 and upper<1
    middle=dist*upper**n
    light=lower**(n-1)/dist
    assert middle<middle_floor and light>light_ceiling
    print("per-step radius, boundary length n, middle/heavy above n excluded, light/middle below n",
          radius,n,middle.upper().str(10),light.lower().str(10))


if __name__=="__main__":
    test(".26",".18",7,arb(".00306")/arb(".8659"),arb(".00000629")/arb(".00306"))
    test(".160",".125",5,arb(".000265")/arb(".01410"),arb(".0000137")/arb(".000265"))
