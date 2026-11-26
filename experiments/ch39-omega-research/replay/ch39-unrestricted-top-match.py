"""Exclude all nonnegative top inputs outside the book band at zero 11431."""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import arb

inverse = run_path(str(Path(__file__).with_name("ch39-top-orbit-inverse.py")))
single = inverse["single"]
source, shell, events, q = single["source"], inverse["shell"], inverse["events"], inverse["q"]
solver = shell["shortest"]
model, grid, ratios = solver["model"], solver["reference"], solver["ratios"]
running, loop = grid["running"], grid["loop"]
anchor = q("26.81666171")
tolerance = arb(2)**-21 + q("1.3e-16")
ports = [[x*x for x in word] for word in source["chord_words"]]
mass = source["weighted_ports"](ports, single["counts"])


def initial(top, upper):
    levels = events["high"] if upper else events["low"]
    state = [arb(x)/mass[0] for x in mass[1:]]+[q(levels[k]) for k in events["names"]]
    state[5] = q(top)**2
    return state


def coefficients(state, order, upper, flow):
    if "NEUTRAL_YUKAWA_BOUNDS" in flow:
        flow = {**flow,"NEUTRAL_YUKAWA":flow["NEUTRAL_YUKAWA_BOUNDS"][int(upper)]}
    return grid["rational"]["taylor"](state, order, flow, 1/loop, opposite_trace=not upper)


def advance(state, step, upper, *, need_tube=False, flow=running, integrate_state=False):
    rate = arb(1)
    for _ in range(20):
        tube = [x+arb(0,(abs(step)*rate).upper()) for x in state]
        derivative = coefficients(tube, 1, upper, flow)[1]
        if all(abs(v) < rate for v in derivative):
            break
        rate *= 2
    else:
        raise ArithmeticError("no invariant step enclosure")
    segment = arb(0).union(step)
    tube = [x+segment*v for x,v in zip(state,derivative,strict=True)]
    terms = coefficients(state,12,upper,flow)
    end = coefficients(tube,13,upper,flow)[-1]
    if integrate_state:
        integral = [v/14 for v in end]
        for k in reversed(range(len(terms))):
            integral = [v/(k+1)+step*w for v,w in zip(terms[k],integral,strict=True)]
        integral = [step*v for v in integral]
    for term in reversed(terms):
        end = [a+step*b for a,b in zip(term,end,strict=True)]
    if integrate_state:
        return end,tube,integral
    return (end,tube) if need_tube else end


@lru_cache(None)
def anchor_state(top, upper):
    state = initial(top,upper)
    for _ in range(192):
        state = advance(state,anchor/192,upper)
    return state


def residuals(lo, hi, time):
    target = model["target"](11431,time)
    pair = []
    for top,upper in ((lo,False),(hi,True)):
        state = advance(anchor_state(top,upper),time-anchor,upper)
        pair.append([1/(r*u)-target for r,u in zip(ratios,state)])
    return pair


def separate(lo, hi):
    time = anchor
    for _ in range(3):
        low, high = residuals(lo,hi,time)
        time = (time-(max(low)+min(high))/(4*model["target"](11431,time))).mid()
    low, high = residuals(lo,hi,time)
    assert 23 < time < 33
    if max(low) > tolerance and min(high) < -tolerance:
        return 1
    assert hi-lo > Fraction(1,100000)
    mid = (lo+hi)/2
    return separate(lo,mid)+separate(mid,hi)


def increasing_gauges():
    # The right sides for X_i=1/g_i^2 are bounded below on X_j >= 1/2.
    # Linear continuation from the initial point closes this bootstrap to 33.
    for i in range(3):
        lower = -2*(running["BETA"][i]+2*sum(max(v,0) for v in running["B"][i])/loop)/loop
        start = arb(mass[0])/mass[i+1]
        assert min(start, start+33*lower) > arb(1)/2
        assert 1+lower > 2*ratios[i]*tolerance
        assert all(v >= 0 for v in running["D"][i])
        assert all(running["B"][i][j] >= 0 for j in range(3) if j != i)
    assert all(v >= 0 for row in running["C"] for v in row)


def high_tail():
    weights = [q(".37"),arb(-1),q(".63")]
    def projected(rows):
        return [sum(w*row[j]/r for w,row,r in zip(weights,rows,ratios,strict=True)) for j in range(3)]
    b,d = projected(running["B"]),projected(running["D"])
    one = sum(w*v/r for w,v,r in zip(weights,running["BETA"],ratios,strict=True))
    norm, loss_integral = arb(0), arb(0)
    state = initial(Fraction(5,4),False)
    step = arb(27)/256
    for i in range(256):
        state, tube = advance(state,step,False,need_tube=True)
        squares = [u**2 for u in tube[:3]]
        norm += step*2/loop**2*max(sum(max(x,0)*u for x,u in zip(row,squares))
                                   for row in running["B"])
        loss_integral += step*sum(max(-x,0)*u for x,u in zip(b,squares))
        if (i+1)*step >= 23:
            traces = [sum(tube[3+j:6+j]) for j in (0,3,6)]
            speed = -2/loop*(one+(sum(x*u for x,u in zip(b,tube))-
                                  sum(x*y for x,y in zip(d,traces)))/loop)
            assert speed < 0
    loss = (4/loop**2)*norm.exp()*loss_integral
    assert min(d) > loss
    boundary = [1/(r*u) for r,u in zip(ratios,state)]
    invariant = sum(w*x for w,x in zip(weights,boundary))
    assert invariant > q(".001") > 2*tolerance
    assert max(boundary)-model["target"](11431,arb(27)) > tolerance
    print("high-top invariant at 27", invariant.lower().str(10),
          "integral loss below",loss.upper().str(10),"minimum matter coefficient",min(d).lower().str(10))


def verify():
    increasing_gauges()
    low, high = (Fraction(".8616")*u for u in (Fraction(".99"),Fraction("1.01")))
    assert separate(Fraction(0),low) == 20
    assert separate(high,Fraction(5,4)) == 20
    high_tail()
    print("all nonnegative initial top couplings outside the original 1% band excluded at 11431")


if __name__ == "__main__":
    verify()
