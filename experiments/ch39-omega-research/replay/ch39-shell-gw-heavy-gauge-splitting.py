"""Flavor-dependent two-loop gauge splitting directly above the first heavy doublet."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

inert = run_path(str(Path(__file__).with_name("ch39-shell-gw-two-inert-heavy.py")))
flow = inert["flow"]
weights = [5,-3,-2]
charges = [[Fraction(17,10),Fraction(3,2),2], [Fraction(1,2),Fraction(3,2),2], [Fraction(3,2),Fraction(1,2),0]]


def check(index):
    heavy, light = inert["certify"](index)
    assert [sum(a*b for a,b in zip(weights,row)) for row in charges] == [0,-6,6]
    interval = (heavy[1]/heavy[0]).log()
    x = interval/flow["loop"]
    start = (heavy[0]/1000).log()/flow["loop"]
    gauge = flow["gauge"](start,(flow["start"],flow["beta"]))
    slopes = [b+d for b,d in zip(flow["beta"],(arb(1)/10,arb(1)/6,0))]
    squared = flow["gauge"](arb(0).union(x),(gauge,slopes))
    # In a two-doublet quark/lepton matrix system the gauge coefficient is scalar.
    assert all(sum(c*v for c,v in zip(row,squared)) < 12 for row in flow["C"])
    cap = arb(4)
    # A single-site sign flip exposes both heavy columns at that site. D0 is
    # conjugate to F0 up to a unit phase, preserving its entire singular spectrum.
    F,D = inert["fit"]["readouts"](index)
    G = inert["fingerprints"](index)[0]
    lower,upper = [],[]
    for matrix in (F,D):
        local = inert["acb_mat"]([row[:9] for row in inert["rows"](matrix)])
        gram = inert["adj"](local)*local
        for limit,k in ((".0374",0),(".52",9)):
            assert inert["count"](inert["rows"](gram),arb(limit)**2)==k
        scale = 2*inert["g"]*G[1,1].real.sqrt()
        lower.append(arb(".0374")*scale)
        upper.append(arb(".52")*scale)
    # All six quark and charged-lepton matrices, on both branches, remain
    # inside the common Frobenius ball by the first-exit bound.
    drift = x*(12*cap+22*cap**3)
    light_norms = [inert["hs"](m,m).real.sqrt() for m in light]
    assert all(v+drift < cap for v in light_norms+upper)
    assert max(light_norms[1:]) < arb(".011") and drift < arb(".012")
    assert lower[1] > drift
    coeff = 6/(32*arb.pi()**3)
    at_start = coeff*lower[1]**2
    assert at_start > arb(".000351")
    # For the finite integral discard positive terms from the inert lepton and
    # active light down sectors. The inert heavy down matrix stays zero.
    excess = (lower[1]-drift)**2-(light_norms[1]+drift)**2-(light_norms[2]+drift)**2-drift**2
    displacement = coeff*interval*excess
    limit = arb("1.9e-7" if index==43 else "3.7e-7")
    assert displacement > limit
    print("zero",index,"log interval",interval.str(9),"initial downward derivative exceeds",
          at_start.lower().str(9),"leading two-loop downward displacement exceeds",displacement.lower().str(9))


if __name__ == "__main__":
    for n in (43,195):
        check(n)
