"""Constrain the conjugate orientation in a two-copy radial Schur response."""
from pathlib import Path
from runpy import run_path
from fractions import Fraction

from flint import arb

model = run_path(str(Path(__file__).with_name("ch39-radial-gain-two-angle-fit.py")))


def verify(response):
    circle = model["circle"]
    shell = circle["shell"]
    inverse = circle["verify"](response)
    gain = arb(18)/25
    circle["exclude_gain_window"](inverse, (-gain).union(arb(7)/10))
    # With total gain 18/25, an odd gain above 7/10 forces the reversed weight below 1/72.
    assert (1-Fraction(7,10)/Fraction(18,25))/2 == Fraction(1,72)
    s,t,u = [v+arb(0,shell["sine_radius"]) for v in shell["centers"][:3]]
    assert s*t-(1-s*s).sqrt()*(1-t*t).sqrt()*u-u > arb(".0019237")

    _, params = model["isolate"](("2.56051","3.0833"), response, gain)
    windows = [shell["sine_radius"]]*3+[shell["phase_radius"]]
    def deviations(r):
        effective = response+r*(response.conjugate()-response)
        w = model["model"]["observables"]((gain,*params),response=effective)
        ud,ub,cd,cb = [w[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2))]
        c13 = (1-abs(ub)**2).sqrt()
        z = -ud*ub.conjugate()*cd.conjugate()*cb-abs(ud*cb*ub)**2/(1-abs(ub)**2)
        values = (abs(w[0,1])/c13, abs(cb)/c13, abs(ub), z.arg())
        return [v-c for v,c in zip(values,shell["centers"])]
    for i in range(64):
        r = (arb(i)/32768).union(arb(i+1)/32768)
        assert all(abs(v)<rad for v,rad in zip(deviations(r),windows))
    assert deviations(arb(1)/256)[3] < -windows[3]
    print("All four windows hold through reversed weight 1/512 at the original fit.")
    print("At 1/256 the CP window fails there; every fit requires reversed weight < 1/72.")


if __name__ == "__main__":
    verify(model["circle"]["two_ticks"]["verify"]())
