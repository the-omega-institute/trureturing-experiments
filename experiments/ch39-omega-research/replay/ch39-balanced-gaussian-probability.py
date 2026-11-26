"""正负行列式在每个共同核尺度与质量处配平时的必要幂数。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

neutral = run_path(str(Path(__file__).with_name("ch39-neutral-clock-response.py")))


def verify():
    ctx.prec, ctx.cap = 768, 3
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    min_x = arb(117387924448)
    second_x = min_x*(first/second)**2
    # 各荷的一、二阶迹相同，三阶以后的对数余项在 [0,r]。
    def log_remainder(x):
        return arb(19)/3*(36/x)**3
    r1, r2 = [log_remainder(x) for x in (min_x, second_x)]
    required = -(arb(".00677850")).log()/(r1+r2)
    assert required > arb("3.20571e27")
    # 第三个节点在算术矩给定的 x 窗内，用全荷 Hadamard 上界除以中性下界。
    x = arb(".00015")
    edges = neutral["shell"]["geo"]["edges"]
    degrees = [sum(i in (a,b) for a,b,_,_ in edges) for i in range(19)]
    delta = arb(23)/(25*2**20)/(first**2+1)
    cap = prod(x+j for j in degrees)/(x*neutral["neutral_polynomial"]()(x)*(1-delta)**3)
    assert cap < 40860 and cap.log() > r1
    theta_max = (second/first).log()/(arb(23672421776)/first).log()
    any_power = arb(259).log()/(r1+r2+theta_max*(cap.log()-r1))
    assert any_power > arb("23.09066")
    print("逆平方先验所需正幂总数下界", required.str(14))
    print("高处的同质量全荷行列式倍率上界", cap.str(14))
    print("任意实幂先验所需正幂总数下界", any_power.str(14))


if __name__ == "__main__":
    verify()
