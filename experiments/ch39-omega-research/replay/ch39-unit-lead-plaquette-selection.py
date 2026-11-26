"""联合改变零点、端口次数和尺度，用原单位边规范端口及有限计次代价作选择。"""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
leads = run_path(str(root/"ch39-unit-lead-gauge-shell.py"))
winding = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))
flow = leads["floating"]


def phase_minimum():
    theta = winding["angle"](11431)
    resonance = 2*arb.pi()/3
    count = 23226
    cost = arb("1e-13").union(arb("1e-12"))
    delta = resonance-count*theta
    assert 0 < delta < theta/2
    # 路径本征值以 g 为单位，消去径向自由度后的能量以 g v^4 为单位。
    energy = 2*(1-delta.cos())
    trial = 3*energy+count*cost
    assert trial < arb("2.53e-8")
    assert 21000*cost > 3*energy
    assert (count+21000)*theta < 2*resonance-delta
    # 先前每步的本征值下降量以末步和起点两个界取小。
    last_drop = 4*(theta/2).sin()*(delta+theta/2).sin()
    early_drop = 4*(theta/2).sin()*resonance.sin()
    # d[3 E/(1+g E/(2 lambda))]/dE，E<=3 且 g/lambda<=1/4。
    slope = 3/(1+arb(3)/8)**2
    assert min(last_drop.lower(), early_drop.lower())*slope > cost.upper()
    print("trial phase and count, normalized upper",trial.upper().str(10),flush=True)
    return trial.upper()


def verify():
    phase = phase_minimum()
    gauge = flow["gauge"]
    time = arb("26.81676714125598")
    eps = arb("1.005e-7")
    assert all(err < eps for err in gauge["gauge_match"](gauge["word"],11431,time,(18,6,3)))
    flow["flow_bounds"]()
    target = gauge["shortest"]["model"]["target"]
    center = target(11431,time)
    d = arb(".00005")
    drift = eps+arb(".09")*d
    near = min((target(11430,time+d)-center-drift).lower(),
               (center-drift-target(11432,time-d)).lower())
    far = arb(".1435")*d-2*eps
    assert near>0 and far>0
    for count in (1,2,3):
        common, contrast = [arb(v) for v in leads["terminal_modes"](count)]
        assert common < contrast
        lower = min((3*common*near**2).lower(), (contrast*far**2/2).lower())
        trial = 3*contrast*eps**2+arb(".0002")*phase
        print("leads",count,"gauge competitor lower",lower.str(10),
              "total trial upper",trial.upper().str(10),"margin",(lower-trial).lower().str(10),flush=True)
        assert lower-trial > arb("4e-12")


if __name__=="__main__":
    verify()
