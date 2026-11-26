"""Companion neutral scale in an SU(2)_X x U(1)_A Higgs completion of the lightest Y=2 vector."""
from pathlib import Path
import runpy
from flint import arb

cross = runpy.run_path(str(Path(__file__).with_name("ch39-spin-gauge-crossing.py")))
q,b,m = cross["q"],cross["b"],cross["m"]
gs = b["model"]["gauge_inputs"]
pi = arb.pi()
for lh,rh,hmax in ((144,4716,30),(216,4068,27)):
    t,z,_,_ = cross["fit"](lh,rh,hmax,gs)
    # Largest threshold logarithm bounds the lightest charged vector. Before it
    # opens, the full-family threshold and all of the right tower are inactive.
    peak = z.upper()/18 + b["shift"].upper()-b["shift"].lower()
    ton = t.lower()-peak
    active_left = m["wl"].lower()-18*peak
    Gpre = (1/gs[0]**2-(m["beta"][0]*ton+q(lh)/90*active_left)/(8*pi**2)).upper()
    gprime = (q("3/5")/Gpre).sqrt().lower()
    # For each nontrivial SU(2) representation used by the vacuum,
    # j(j+1)-m^2 >= j >= 1/2. Hypercharge neutrality imposes X=-2m.
    ratio = ((4*pi).sqrt()/gprime).upper()
    cutoff = ratio*(-z.lower()/18).exp()
    assert cutoff<1
    print(lh,"gY at first right threshold >",gprime.str(12),"mass ratio <",ratio.str(12),
          "neutral mass / Q <",cutoff.upper().str(12))
