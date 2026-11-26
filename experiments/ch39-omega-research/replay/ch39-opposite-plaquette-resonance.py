"""负取向条件平衡的共振、两种手征的热赔率与两条边的反射。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
endpoints = run_path(str(root/"ch39-classical-endpoint-ensemble.py"))
q = arb


def oscillator_free_over_heat(count, temperature_ratio):
    unit = q.pi()/(2*count)
    return 6*sum((2*((j*unit).sin()/temperature_ratio).sinh()).log() for j in range(1,count))


def verify():
    ctx.prec = 128
    spread = endpoints["volume_bounds"]()
    positive, negative, tail_start = 23226, 46453, 22000
    heat_min, heat_max = q("1.25e-18"), q("2e-17")
    phase = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))
    theta = phase["angle"](11431)
    alpha = 2*q.pi()/3
    delta = 2*alpha-negative*theta
    chord = lambda z: 2*(1-z.cos())
    e1 = chord(alpha-positive*theta)/positive
    e2 = chord(delta)/negative
    assert -theta/2 < delta < 0 and e2>e1
    assert (negative+tail_start)*theta < 2*q.pi()
    slope = 3/(1+q(3)/8)**2
    early = slope*4*(theta/2).sin()*(delta+theta/2).sin()/negative-q("2e-16")
    middle = slope*(chord(theta-delta)/(negative+tail_start)-e2)
    late = tail_start*q("1e-17")-3*e2
    base_exponent = min(v.lower() for v in (early,middle,late))/heat_max
    assert (-q(".5")).exp() < q(2)/3 and early/heat_max > 1
    lengths = (tail_start+5)*(-base_exponent+spread).exp()
    # 第一段相位位于上半圆，负取向的受限势至少为 2/p。
    oscillator_top = 12*negative*q("4e-17")/q.pi()
    wrong_arc = 16*negative*(spread+(-q(1)/negative+3*e2+oscillator_top)/heat_max).exp()

    equilibrium = run_path(str(root/"ch39-thermal-plaquette-occupation.py"))
    joint = equilibrium["joint"]
    target = joint["flow"]["gauge"]["shortest"]["model"]["target"]
    time, probe, tube = q("26.81676714125598"),q("1e-8"),q(".00005")
    center,eps = target(11431,time),q("1.005e-7")
    allowance = eps+(q(".09")+2*center*(2*probe).exp())*probe
    split_scale = q(".1435")*tube-2*eps
    separation = min((target(11430,time+tube)-center).lower(),
                     (center-target(11432,time-tube)).lower())-eps-q(".09")*tube
    assert separation>0 and split_scale>0
    other_labels = []
    witnesses = []
    for leads in (1,2,3):
        common,contrast = (q(v) for v in joint["leads"]["terminal_modes"](leads))
        lower = min((3*common*separation**2).lower(),(contrast*split_scale**2/2).lower())
        witness = 3*e2+oscillator_top+3*contrast*allowance**2
        exponent = (lower-witness)/heat_max-spread-16-(15/probe).log()
        other_labels.append(exponent.lower())
        witnesses.append(witness.upper())
    error = lengths+wrong_arc+(-min(other_labels)).exp()
    assert -error.log()>2160
    print("负取向编号或链长错配的概率负对数下界",(-error.log()).lower().str(10))

    # 取向翻转前后共享的两条未改路径，若两端都规则，其有向相位符号一致。
    barrier = 4*q(".02")**2
    bottleneck = (barrier/2-max(witnesses))/heat_max-spread-16-(16*15/probe).log()
    assert bottleneck>q("3.99e13")

    # 两个主扇区在同一编号，尺度积分和编号先验相消。
    assert negative == 2*positive+1
    lower_free = oscillator_free_over_heat(negative,q(1)/2)-oscillator_free_over_heat(positive,q(1)/2)
    upper_free = oscillator_free_over_heat(negative,q(1)/4)-oscillator_free_over_heat(positive,q(1)/4)
    gap_min = 3*(e2-e1)/((1+e1/8)*(1+e2/8))
    gap_max = 3*(e2-e1)/((1+e1/12)*(1+e2/12))
    # 两侧条件集中误差均 < 1/2，主扇区权重与全取向权重的比值在 [1/2,1]。
    log_low = lower_free+gap_min/heat_max-spread-q(2).log()
    log_high = upper_free+gap_max/heat_min+spread+q(2).log()
    print("正/负取向概率比的自然对数窗",log_low.lower().str(10),log_high.upper().str(10))
    assert log_low > 139000 and log_high < 341000

    # 翻转两条边的相位并将奇长度减半。原负取向基模上的路径系数下降。
    split2 = chord(negative*theta-alpha)/negative-e2
    assert e1 < e2 and 4/positive < 3*split2 and 3*(6-e2)>4/positive
    assert lower_free/6-spread-q(16).log()>23170
    assert error+(spread+q(16).log()-lower_free/6).exp() < (-q(2160)).exp()


if __name__ == "__main__":
    verify()
