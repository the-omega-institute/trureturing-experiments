"""Evaluate normalized gauge residuals at individual finite input and clock codes."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

volume = run_path(str(Path(__file__).with_name("ch39-volume-contenders.py")))
model, q, running = volume["gauge"], volume["q"], volume["running"]
loop = model["loop"]


rational = run_path(str(Path(__file__).with_name("ch39-rational-gauge-flow.py")))


def taylor(initial, order):
    return rational["taylor"](initial, order, running, 1/loop)


def advance(initial, step, degree):
    # An invariant enclosure for the portion of the solution on this step.
    pad = abs(step)/2
    tube = [x+arb(0,pad) for x in initial]
    assert all(abs(f)*abs(step) < pad for f in taylor(tube,1)[1])
    endpoint = taylor(initial,degree-1)
    remainder = taylor(tube,degree)[-1]
    for term in reversed(endpoint):
        remainder = [a+step*b for a,b in zip(term,remainder)]
    return remainder


def readout(codes, clock, yukawa):
    assert len(codes) == 3 and all(type(k) is int and 0 <= k < 2**32 for k in codes)
    assert type(clock) is int and 0 <= clock < 2**45
    assert len(yukawa) == 9
    low, high, _ = volume["input_bounds"]()
    values = [1/(a+(b-a)*(2*k+1)/2**33) for a,b,k in zip(low,high,codes)]
    values += [y**2 for y in yukawa]
    time = 23 + arb(10)*(2*clock+1)/2**46
    steps = 256
    for _ in range(steps):
        values = advance(values,time/steps,14)
    return time,[1/(u*r) for u,r in zip(values,model["ratios"])]


def verify():
    ctx.prec = 160
    codes = [382269474,3503200819,2112544438]
    clock = 13349175109351
    time, end = readout(codes,clock,model["yukawa_centers"])
    for label in (11736,11737,11738):
        target = model["target"](label,time)
        residuals = [v-target for v in end]
        assert max(v.rad() for v in residuals) < arb(2)**-42
        worst = max(abs(v).upper() for v in residuals)
        if label == 11737:
            assert worst < arb(2)**-36
        else:
            assert min(abs(v).lower() for v in residuals) > arb(2)**-14
        print(label,"normalized residuals",*[v.str(7) for v in residuals],
              "maximum radius",max(v.rad() for v in residuals).str(4))


if __name__ == "__main__":
    verify()
