"""A comparison bound for the fully running quark matrices with tree shell matching at each mass."""
from pathlib import Path
from runpy import run_path

from flint import arb

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
fit = flow["fit"]


def check(index):
    Q,heavy = flow["masses"](index)
    x1 = (heavy[0]/1000).log()/flow["loop"]
    xQ = (Q/1000).log()/flow["loop"]
    gfirst = flow["gauge"](x1,(flow["start"],flow["beta"]))
    unified = flow["gauge"](xQ,(flow["start"],flow["beta"]))[2]
    # The inverse squared gauge couplings are concave in log(mu), since each
    # new doublet raises or preserves b. Minima on [first,Q] occur at an end.
    caps = [arb.max(v,unified) for v in gfirst]
    damping = arb(9)
    assert all(2*sum(a*v for a,v in zip(c,caps)) < damping for c in flow["C"][:2])
    _,_,_,_,F,D,lower,*_ = flow["at_threshold"](index)
    strength2 = lower.lower()**2
    total = strength2*(fit["hs"](F,F).real+fit["hs"](D,D).real)
    tlast = (heavy[-1]/heavy[0]).log()/flow["loop"]
    S = total*(-damping*tlast).exp()
    assert S>damping
    remaining = flow["loop"]/damping*(S/(S-damping)).log()
    pole = heavy[-1]*remaining.exp()
    assert pole < Q
    # In either 9-dimensional flavor block a trace forces a large partial wave.
    cap = 48*arb.pi()
    assert S < cap
    delay = flow["loop"]/damping*((1-damping/cap)/(1-damping/S)).log()
    wave_scale = heavy[-1]*delay.exp()
    assert wave_scale < pole
    print("zero",index,"last threshold",heavy[-1].str(8),"quark squared norm above",S.lower().str(8),
          "running contact-limit partial-wave criterion crossed by",wave_scale.upper().str(10),
          "one-loop blow-up by",pole.upper().str(10),"gauge-only crossing",Q.str(9))


if __name__ == "__main__":
    for n in (43,195):
        check(n)
