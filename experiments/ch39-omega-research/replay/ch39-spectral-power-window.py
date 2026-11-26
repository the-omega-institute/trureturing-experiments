"""四个时间目标的有限高度窗对幂律谱权的限制。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


shared = run_path(str(Path(__file__).with_name("ch39-shared-flavor-clock-weights.py")))


def verify():
    ctx.prec, ctx.cap = 640, 3
    times = shared["clock"]["clock"]["times"]()
    moments = [v/200**i for i, v in enumerate(shared["clock"]["clock"]["cleared_moments"](times, 4))]
    a, b = arb(".00015"), arb(".015")
    residual = moments[1] - a*moments[0] - moments[4]/b**3
    assert residual > 0
    mass = residual**2/moments[2]
    assert mass > arb(".29045")
    first = acb.zeta_zero(8359).imag
    cap = arb("2.634163e-9")
    ratio = (b/a).sqrt()
    assert (arb(41312660)/b.sqrt()) > first and first > 8000
    # 对 H/R≤γ≤H，Σ1/γ ≤ (1+log R)log H/(2π)。
    coefficient = cap*first*(1+ratio.log())/(2*arb.pi())
    min_log_height = mass/coefficient
    min_log_scale = min_log_height + a.log()/2
    assert min_log_scale/arb(10).log() > 10764
    # u exp(-δu) ≤ 1/(eδ)。此式只在 δ>0 时使用。
    max_delta = 1/(arb(1).exp()*(mass/coefficient-(first*ratio).log()))
    assert max_delta < arb(".00001485")
    # 逆高度平方的整带总量。标准计数主项求 Stieltjes 积分，两个端点保留误差。
    low = arb("41312660.31094738")/b.sqrt()
    high = low*ratio
    def main_tail(t):
        return ((t/(2*arb.pi())).log()+1)/(2*arb.pi()*t)
    count_error = arb(".28")*((2*low.log()+arb(".5"))/low**2 + high.log()/high**2)
    upper_band = main_tail(low) - main_tail(high) + count_error
    gain = mass/(cap*first**2*upper_band)
    assert gain > 195557000
    print("相对逆平方谱，十倍高度带内的最低效率倍率", gain.str(14))
    print("所需有限窗权重、log10(A) 下界、幂指数超过一的最大余量",
          mass.str(12), (min_log_scale/arb(10).log()).str(12), max_delta.str(12))


if __name__ == "__main__":
    verify()
