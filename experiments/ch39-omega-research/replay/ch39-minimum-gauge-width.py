"""Minimum displacement of the low-energy gauge inputs needed for exact matching."""
from pathlib import Path
from flint import arb, ctx
from runpy import run_path

model = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))
q = model["q"]


def obstruction():
    delta = q(".000618")
    scale = (1-delta).union(1+delta)
    inverse = [1/(scale*v)**2 for v in model["centers"]]
    logs = [2*(v*q(".99").union(q("1.01"))).log() for v in model["yukawa_centers"]]
    time = q("26.81470")
    end,_,weak_minus_color = model["evolve"](inverse,logs,time,16384)
    hyper_minus_weak = end[0]/model["ratios"][0]-end[1]/model["ratios"][1]
    assert hyper_minus_weak>q(".000004") and weak_minus_color<q("-.000009")
    def slopes(step, state, rate):
        if step>=23*64:
            normalized = [v/r for v,r in zip(rate,model["ratios"])]
            assert normalized[0] < normalized[1] < normalized[2]
    model["running"]["integrate"]([scale*c for c in model["centers"]],33*64,slopes)
    print("618 ppm exclusion, two separating margins",hyper_minus_weak.lower().str(8),
          (-weak_minus_color).lower().str(8))


def first_witness():
    n = 11462
    time = q("26.81444904")
    common = model["target"](n,time)
    logs = list(map(q,("-25.412925110863 -12.961070951161 -1.515705670913 "
                       "-23.787163242603 -17.810570965491 -10.094520128717 "
                       "-25.620344619829 -14.919786191565 -9.254900631546").split()))
    inv,y,_ = model["evolve"]([common*r for r in model["ratios"]],logs,-time,131072,
                               diagonal_quarks=True)
    shifts = [(1/(v.sqrt()*c)-1)*10**6 for v,c in zip(inv,model["centers"])]
    assert max(abs(x).upper() for x in shifts)<631
    assert all(q(".99")<(v/2).exp()/c<q("1.01") for v,c in zip(y,model["yukawa_centers"]))
    print("631 ppm witness with zero",n,"relative input shifts in ppm",
          *[(float(v.lower()),float(v.upper())) for v in shifts])


if __name__ == "__main__":
    ctx.prec = 144
    obstruction()
    first_witness()
