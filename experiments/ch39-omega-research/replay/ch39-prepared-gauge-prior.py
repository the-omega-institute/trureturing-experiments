"""An everywhere-positive smooth gauge preparation reversing the 11737 volume ranking."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
volume = run_path(str(root / "ch39-volume-contenders.py"))
gauge, q = volume["gauge"], volume["q"]
ratios = gauge["ratios"]


def verify():
    ctx.prec = 112
    fixed_y = [2*(c*q(".993")).log() for c in gauge["yukawa_centers"]]
    label = 10875
    time = q("26.859482438")
    value = gauge["target"](label,time)
    # The old six-face witness has a uniform signed residual margin > 2e-6
    # and the maximum row sum of the endpoint gauge Jacobian is < 1.023456.
    radius = q("1e-7")
    movement = q("1.023456")*radius
    assert all(movement/r < q("2e-6") for r in ratios)
    assert all(q(".99") < (y/2).exp()/c < q("1.01") for y,c in zip(fixed_y,gauge["yukawa_centers"]))

    factors = q(".999").union(q("1.001"))
    def observe(i, inverse, rate):
        if i >= 23*64:
            velocities = [v/r for v,r in zip(rate,ratios)]
            assert velocities[2]-velocities[0] > q(".12")
            assert all(abs(v) < q(".16") for v in velocities)
    gauge["running"]["integrate"]([c*factors for c in gauge["centers"]],33*64,observe)
    tolerance = arb(2)**-21
    displacement = (2*tolerance + movement*(1+1/ratios[0]))/q(".12")
    assert displacement < q(".00001")
    # A match must be close to the one time at which the reference input has
    # equal hypercharge and color. Every other target is outside the color band.
    span = time + arb(0,displacement)
    color = value + arb(0,movement+q(".16")*displacement+tolerance)
    assert gauge["target"](label-1,span) > color
    assert gauge["target"](label+1,span) < color
    half_time = arb(2)**-25
    local_target = gauge["target"](label,time+arb(0,half_time))
    assert all(movement/r + (2*local_target+q(".16"))*half_time < tolerance for r in ratios)

    targets, _ = volume["zero_targets"](tolerance)
    error = q(".000097")/10000**2
    min_weight = targets[label][1]-error
    max_weight = max(w.upper() for _,w in targets.values())+error
    assert 0 < max_weight < 2*min_weight
    # The narrow component is any unit-mass smooth bump supported in the
    # certified cube. Add uniform background to make its density positive.
    background = arb(2)**-30
    odds = (1-background)*(2*half_time/10)*min_weight/(background*max_weight)
    assert odds > q("6.3893")
    print("10875 isolated in gauge cube of radius",radius.str(4),"matching time displacement <",
          displacement.upper().str(9),"relative score for either protocol >",odds.lower().str(9))


if __name__ == "__main__":
    verify()
