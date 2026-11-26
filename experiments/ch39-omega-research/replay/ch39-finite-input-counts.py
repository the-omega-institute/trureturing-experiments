"""Transfer the shell-weighted volume certificate to a finite equal-count input grid."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    ctx.prec = 112
    volume = run_path(str(Path(__file__).with_name("ch39-volume-contenders.py")))
    model, q = volume["gauge"], volume["q"]
    ratios = model["ratios"]
    half_width = arb(2)**-22
    broadest = arb(2)**-21
    gauge_bits, time_bits, readout_bits = 32, 45, 40
    x_low,x_high,_ = volume["input_bounds"]()
    max_width = max(((b-a)/min(ratios)).upper() for a,b in zip(x_low,x_high))
    # At fixed initial Y the previously certified row norm of the gauge flow is
    # at most 1.023456 on the convex low-input cube, for all relevant times.
    displacement = q("1.023456")*max_width/2**(gauge_bits+1)
    factors = q(".999").union(q("1.001"))
    def check(i, inverse, rate):
        if i >= 23*64:
            for x, v, r in zip(inverse,rate,ratios):
                assert 0 < x/r < 6
                assert abs(v/r) < q(".16")
    model["running"]["integrate"]([c*factors for c in model["centers"]],33*64,check)
    clock_step = arb(10)/2**(time_bits+1)
    # At either matching endpoint G_n < 6+eps+roundoff. Propagate its exact
    # exponential time dependence over the cell, independently of the zero label.
    clock_error = (6+broadest)*((2*clock_step).exp()-1) + q(".16")*clock_step
    total = displacement + clock_error + arb(2)**(-readout_bits)
    assert 0 < half_width-total < half_width+total < broadest
    comparison = []
    for dim, gap in ((2, ".0001176"), (3, ".00007643")):
        advantage = (1+q(".0002"))*((half_width-total)/(half_width+total))**dim-1
        assert advantage > q(gap)
        comparison.append(advantage.lower())
    print("Input bits",gauge_bits,"per gauge; clock bits",time_bits,"residual error at most 2^-",readout_bits)
    print("Cell and readout uncertainty <",total.upper().str(9),"finite-count relative advantages >",
          *[x.str(9) for x in comparison])


if __name__ == "__main__":
    verify()
