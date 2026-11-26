"""将端口的 p 次重复展开为单位相位链，检验 Schur 刚度和联合最小值。"""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
joint = run_path(str(root/"ch39-unit-lead-plaquette-selection.py"))
q = arb


def verify():
    count = 23226
    theta = joint["winding"]["angle"](11431)
    resonance = 2*q.pi()/3
    delta = resonance-count*theta
    cost = q("1e-17").union(q("2e-16"))
    assert 0<delta<theta/2
    energy = lambda angle: 2*(1-angle.cos())
    phase = 3*energy(delta)/count
    trial = phase+count*cost
    assert 21000*cost>phase
    # 区间内到最近的非平凡三次单位根至少有以下角距。
    next_distance = min((theta-delta).lower(),(2*resonance-(count+21000)*theta).lower())
    assert next_distance > delta and energy(next_distance)/(count+21000) > energy(delta)/count
    last_step = 4*(theta/2).sin()*min((delta+theta/2).sin().lower(),resonance.sin().lower())
    radial_slope = 3/(1+q(3)/8)**2
    assert radial_slope*last_step/count > cost.upper()

    eps = q("1.005e-7")
    separation = q(".00005")
    flow = joint["flow"]
    gauge = flow["gauge"]
    time = q("26.81676714125598")
    target = gauge["shortest"]["model"]["target"]
    center = target(11431,time)
    outside = min((target(11430,time+separation)-center).lower(),
                  (center-target(11432,time-separation)).lower())-eps-q(".09")*separation
    spread = q(".1435")*separation-2*eps
    assert outside>0 and spread>0
    for lead_count in (1,2,3):
        m,d = [q(v) for v in joint["leads"]["terminal_modes"](lead_count)]
        lower = min((3*m*outside**2).lower(),(d*spread**2/2).lower())
        assert lower-trial-3*d*eps**2>q("4.4e-12")
    mass = 2*(q.pi()/(2*count)).sin()
    print("normalized phase upper",phase.str(10),"trial upper",trial.upper().str(10),
          "smallest Dirichlet mass",mass.str(12),"unit for twice triplet",(2/mass).str(12),flush=True)


if __name__=="__main__":
    verify()
