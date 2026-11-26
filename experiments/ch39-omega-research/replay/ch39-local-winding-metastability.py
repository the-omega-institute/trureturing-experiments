"""有限 Gibbs 扇区的局部计次动力学：跨共振瓶颈与从短链出发的制备。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
osc = run_path(str(root/"ch39-winding-chain-vacuum-energy.py"))
joint = osc["chain"]["joint"]
q = arb


def verify():
    ctx.prec = 144
    first, second = 23226,46453
    boundary = 34839
    theta = osc["angle"](11431)
    radial=(q(1)/6).union(q(1)/4)
    max_heat = q("2e-17")
    gauge = joint["flow"]["gauge"]
    target = gauge["shortest"]["model"]["target"]
    time, width = q("26.81676714125598"), q("1e-8")
    center,eps = target(11431,time), q("1.005e-7")
    inner=eps+(q(".09")+2*center*(2*width).exp())*width
    spread=q(".1435")*q(".00005")-2*eps
    separation=min((target(11430,time+q(".00005"))-center).lower(),
                   (center-target(11432,time-q(".00005"))).lower())-eps-q(".09")*q(".00005")
    high=osc["phase_energy"](boundary+1,theta,radial)
    assert high>q("8e-5")
    start=osc["phase_energy"](second,theta,radial)+q("4e-17")*12*second/q.pi()
    height = acb.zeta_zero(11431).imag
    total_weights=(1+q.const_euler()/2-(4*q.pi()).log()/2)*(q(1)/4+height**2)
    for leads in (1,2,3):
        common, contrast=[q(v) for v in joint["leads"]["terminal_modes"](leads)]
        cost_other=min((3*common*separation**2).lower(),(contrast*spread**2/2).lower())
        margin=min(high.lower(),cost_other)-start-3*contrast*inner**2
        exponent=margin/max_heat-(total_weights*5/width).log()
        print("leads",leads,"cut versus second basin: log weight ratio upper",(-exponent).str(12),flush=True)
        assert exponent>95000
    # n 固定时的向下能差直接控制 Metropolis 的反向接受率。
    delta=2*q.pi()/3-first*theta
    drop=4*(theta/2).sin()*min((delta+theta/2).sin().lower(),(2*q.pi()/3).sin().lower())
    one_step=3*drop/(1+q(3)/8)**2/first-q("2e-16")
    assert one_step/max_heat>40000
    # 4p 次公平方向提议中，正向不足 p 次的 Hoeffding 上界。
    proposals=4*first
    assert q(proposals)/8>11000
    print("short-chain forward energy gain in thermal units", (one_step/max_heat).str(12),flush=True)


if __name__=="__main__":
    verify()
