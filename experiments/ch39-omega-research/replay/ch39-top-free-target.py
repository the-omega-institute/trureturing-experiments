"""Infer the top input before choosing any common gauge target or zero label."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

comp = run_path(str(Path(__file__).with_name("ch39-unrestricted-top-match.py")))
q, ratios, running, loop = (comp[k] for k in ("q","ratios","running","loop"))
anchor, tolerance = comp["anchor"], comp["tolerance"]


def normalized(top, upper, time):
    state = comp["advance"](comp["anchor_state"](top,upper), time-anchor, upper)
    return [1/(r*u) for r,u in zip(ratios,state)]


def separate(lo, hi):
    time = anchor
    for _ in range(3):
        low = normalized(lo, False, time)
        high = normalized(hi, True, time)
        mid = [(x+y)/2 for x,y in zip(low,high,strict=True)]
        time = (time+(mid[0]-mid[2])/q(".1436")).mid()
    assert 23 < time < 33
    low, high = normalized(lo,False,time), normalized(hi,True,time)
    above = min(low[0]-high[1], low[2]-high[1])
    below = min(low[1]-high[0], low[1]-high[2])
    if max(above,below) > 2*tolerance:
        return 1
    assert hi-lo > Fraction(1,1000000)
    middle = (lo+hi)/2
    return separate(lo,middle)+separate(middle,hi)


def projections(first):
    weights = [q(first),arb(-1),1-q(first)]
    def project(rows):
        return [sum(w*row[k]/r for w,r,row in zip(weights,ratios,rows,strict=True)) for k in range(3)]
    one = sum(w*b/r for w,r,b in zip(weights,ratios,running["BETA"],strict=True))
    return weights, project(running["B"]), project(running["D"]), one


def check_envelopes():
    # The lower envelope at y_t=1.1 controls the projected high tail. The
    # upper envelope at the same endpoint bounds all smaller matter traces.
    high = comp["initial"](Fraction("1.1"), True)
    low = comp["initial"](Fraction("1.1"), False)
    step, cut = arb(1)/16, q("429/16")
    tests = [projections(".37"),projections(".35")]
    loss_integrals = [arb(0),arb(0)]
    norm = arb(0)
    for i in range(33*16):
        start = low
        low, tube = comp["advance"](low,step,False,need_tube=True)
        high, upper_tube = comp["advance"](high,step,True,need_tube=True)
        traces = [sum(tube[k:k+3]) for k in (3,6,9)]
        squares = [u**2 for u in tube[:3]]
        norm += 2*step/loop**2*max(sum(max(x,0)*u for x,u in zip(row,squares))
                                   for row in running["B"])
        for idx,(weights,b,d,one) in enumerate(tests):
            if idx == 0 and i*step >= cut:
                continue
            loss_integrals[idx] += step*sum(max(-x,0)*u for x,u in zip(b,squares))
            if i*step < 23 or (idx == 1 and i*step < cut):
                continue
            speed = -2/loop*(one+(sum(x*u for x,u in zip(b,tube))-sum(x*y for x,y in zip(d,traces)))/loop)
            initial = sum(w/(r*u) for w,r,u in zip(weights,ratios,start))
            assert initial+arb(0).union(step)*speed > q(".0005") > 2*tolerance
        if i*step >= 23:
            # All inverse gauge coordinates exceed 1/2 by the bootstrap,
            # and comparison bounds all nine nonnegative squared Yukawas.
            gauge = [arb(0).union(arb(2))]*3
            matter = [arb(0).union(sum(upper_tube[k:k+3]).upper()) for k in (3,6,9)]
            velocities = [-2/(loop*r)*(b+(sum(a*g-c*y for a,g,c,y in zip(br,gauge,dr,matter,strict=True)))/loop)
                          for r,b,br,dr in zip(ratios,running["BETA"],running["B"],running["D"],strict=True)]
            assert velocities[0] < velocities[1] < velocities[2]
    for (_,_,d,_),integral in zip(tests,loss_integrals,strict=True):
        assert 4*norm.exp()*integral/loop**2 < min(d)


def verify():
    comp["increasing_gauges"]()
    check_envelopes()
    left = separate(Fraction(0),Fraction(".8562"))
    right = separate(Fraction(".8581"),Fraction("1.1"))
    print("arbitrary common target forces 0.8562 < y_t < 0.8581; covering intervals",left,right)


if __name__ == "__main__":
    verify()
