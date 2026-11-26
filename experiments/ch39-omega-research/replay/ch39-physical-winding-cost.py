"""将三态费米子真空能、固定端点弹性能和三条振子链放在同一体积中。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

free = run_path(str(Path(__file__).with_name("ch39-free-winding-cost.py")))


def verify():
    ctx.prec = 160
    lower, upper = free["cost_window"]()
    theta = free["absolute"]["angle"](11431)
    alpha = 2*arb.pi()/3
    p, rival = 23226,92905
    delta = alpha-p*theta
    scalar = 6*(1-delta.cos())/p
    rival_angle = abs(rival*theta-4*alpha)
    assert 6*(1-rival_angle.cos())/rival < scalar
    # h = hbar*Omega / (cell_volume*g*f_e^4)，u = hbar*Omega / (cell_volume*M^4)。
    largest_g = 4*arb.pi()*arb("5e-10")**2
    small_u = arb("1e-17")*largest_g*(arb(1000)/42)**4
    largest_step = 3*arb.pi()/2
    assert small_u*largest_step < lower
    minimum_edge_scale = 42*(lower/(largest_step*arb("1e-17")*largest_g))**(arb(1)/4)

    relative_unit = arb("1e-18")*(arb(10000)/42)**4
    u = arb("1e-17")*relative_unit
    assert arb("1e-18") < largest_g
    smallest_step = 3*arb(2).log()
    heat_upper = u/2
    early_gap = upper-u*largest_step
    # p+1 至 p+20999 的最近角距在两个端点处取得最小。
    cutoff_distance = min((theta-delta).lower(), (2*alpha-(p+21000)*theta).lower())
    assert cutoff_distance > delta
    assert 6*(1-cutoff_distance.cos())/(p+21000) > scalar
    middle_gap = free["vacuum"](cutoff_distance)-free["vacuum"](delta)
    late_slope = u*smallest_step-lower-relative_unit*scalar/21000
    assert early_gap>0 and middle_gap>0 and late_slope>0
    opposite_gap = free["vacuum"](abs((p*theta+alpha).cos().acos()))-free["vacuum"](delta)
    # 每个非目标 p 最多两种取向。晚段从第 21000 项起求几何和。
    a = (-early_gap/heat_upper).exp()
    b = (-late_slope/heat_upper).exp()
    failure = ((-opposite_gap/heat_upper).exp() + 2*a/(1-a)
               + 2*20999*(-middle_gap/heat_upper).exp() + 2*b**21000/(1-b))
    assert failure.log() < -41000
    print("小标尺最大每步费用", (small_u*largest_step).str(11),
          "保住目标所需 f_e/Q_* 下界",minimum_edge_scale.str(12),
          "大标尺 hbar*Omega/(V M^4)",u.str(11),
          "晚段每步能隙下界",late_slope.str(11),"错配概率对数上界",failure.log().str(12))


if __name__ == "__main__":
    verify()
