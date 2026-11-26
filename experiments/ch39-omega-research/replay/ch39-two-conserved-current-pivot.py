"""只读守恒的左弱与颜色电流，对任意允许的矩阵背景筛选全部零点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

absolute = run_path(str(Path(__file__).with_name("ch39-fermion-absolute-matching.py")))


def verify(*, label=11431, center="26.81676714125598", word=None, discarded=(18,6,3), evolved_state=None):
    ctx.prec = 160
    gauge = absolute["gauge"]
    absolute["ground"]["flow_bounds"]()
    model = gauge["shortest"]["model"]
    ratios = gauge["shortest"]["ratios"]
    time = arb(center)
    mass_ratio = 42*(arb("26.81676714125598")-time).exp()
    source = gauge["word"] if word is None else word
    target = model["target"](label,time)
    state = evolved_state if evolved_state is not None else gauge["gauge_state"](source,time,discarded)
    mixing = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    residual = [1/(v*y)-target+arb(0,m.upper()) for y,v,m in zip(state,ratios,mixing)]
    tol = arb(2)**-21
    budget = abs(residual[1]-residual[2]).upper()+2*tol
    speed = arb(".089")-arb(".038")
    prefactor = 1/(4*arb.pi()**2)
    r = arb("4e-5")
    low_m,high_m = 1-3*r,1+3*r
    color_mass = (1-r).union(1+r)
    weak_mass = low_m.union(high_m)
    nominal = (arb(23)-time).union(arb(33)-time)
    scale = nominal
    # K(z) <= log(1+z/4) and K(z)<=z/5. Both subtracted thresholds are positive.
    widths = []
    for _ in range(3):
        z = (scale.upper().exp()/(mass_ratio*low_m))**2
        raw = prefactor*min((1+z/4).log().upper(),(z/5).upper())
        width = (budget+raw*max(1/ratios[1],1/ratios[2]))/speed
        scale = arb(0,width.upper())
        widths.append(width.upper())
    def kernel(z):
        return (z/5-3*z**2/140).lower().union((z/5).upper())
    def subtracted(mass, interval):
        return prefactor*(kernel((interval.exp()/(mass_ratio*mass))**2)
                          -kernel(((-time).exp()/(mass_ratio*mass))**2))
    weak = subtracted(weak_mass,scale)
    color = subtracted(color_mass,scale)
    radius = (budget+abs(weak/ratios[1]-color/ratios[2]).upper())/speed
    scale = arb(0,radius.upper())
    measured = target+residual[2]+arb(".089").union(arb(".09"))*scale-subtracted(color_mass,scale)/ratios[2]
    gap_low = model["target"](label-1,time+radius)-measured-tol
    gap_high = measured-tol-model["target"](label+1,time-radius)
    assert gap_low>0 and gap_high>0
    print("Weak/color-only scale radius",radius.upper().str(10),
          "adjacent target margins",gap_low.lower().str(10),gap_high.lower().str(10))


if __name__ == "__main__":
    verify()
