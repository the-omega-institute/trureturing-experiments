"""端点和链内场同时量子化时的有限链热迹与环取向。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
q = arb


def packet_cost(k_upper, inertia_ratio):
    # d=eta h²/4，h=ℏΩ/(gv⁴)，高斯宽度 u=h sqrt(eta)/(2 sqrt(14+kmax))。
    hmax = q("4e-17")
    heat_max = q("2e-17")
    scale = (14+k_upper).sqrt()
    width_per_h = inertia_ratio.sqrt()/(2*scale)
    # 小盒的每条边模长平方不超过 1+0.09 sqrt(H)。
    linear = 6*(14+k_upper)+6*(26+k_upper)*q(".09")*heat_max.sqrt()
    quadratic = 78+3*k_upper
    return 4*(6*inertia_ratio/(4*width_per_h)+linear*width_per_h
              +quadratic*hmax*width_per_h**2)


def verify():
    ctx.prec = 128
    volume = run_path(str(root/"ch39-classical-endpoint-ensemble.py"))["volume_bounds"]()
    # 下面固定 eta=1，即按无量纲场坐标给两种场同一动能。
    phase = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))
    p = 23226
    phi, alpha = p*phase["angle"](11431), 2*q.pi()/3
    e = 2*(1-(phi-alpha).cos())/p
    reverse = 2*(1-(phi+alpha).cos())/p
    slope_chi = 3*(reverse-e)/(2*(1+e/8)*(1+reverse/8))
    assert slope_chi > 0
    for k_upper, probability_power in ((2,3883),(100,3722),(10000,1577)):
        shift = volume+packet_cost(q(k_upper),q(1))
        assert 4000-shift > probability_power
        deviation = 2*q("2e-17")*(shift+q(16).log())/slope_chi
        print("k 上端",k_upper,"端点热迹对数跨度",shift.str(11),
              "正取向均值缺额上界",deviation.str(11))
        if k_upper==2:
            assert deviation < q("2.48e-11")
    # 质量较重时，在原最窄的 64H 判向窗内也保有足够的自由能差。
    small = volume+packet_cost(q(2),q(".01"))
    assert small<31
    budgets = [(q(200000)/8,10001),(q(800000)/8,50001),(q(1000000)/8,24385),
               (q(299000000)*q("2e-8")/(8*q(".0005")),661)]
    risk = sum((-(mu-target)**2/(2*mu)).exp() for mu,target in budgets)
    risk += 300000000*((small-64).exp()+(-q(190)).exp()+(small-12000).exp())
    assert risk<q("1e-6")


if __name__ == "__main__":
    verify()
