"""用有限振子配分函数检验热平衡的链长选择及高温加长反例。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

osc = run_path(str(Path(__file__).with_name("ch39-winding-chain-vacuum-energy.py")))
q = arb


def verify():
    ctx.prec = 144
    # t log(1/t) 是凹函数，闭区间的下界由端点给出。
    low,high = q(1)/4,q(1)/2
    lower = 6*min((-low*low.log()).lower(),(-high*high.log()).lower())
    upper = 3*q.pi()/2
    assert lower*q("5e-18") > q("1e-17")
    assert upper*q("4e-17") < q("2e-16")

    # 在两倍链长处偶数频率重合，奇数频率用有限正弦乘积估计。
    # 取四倍长度作实际反例，避免依靠无上限的链长。
    p = 23226
    t = q("1.1")
    minus = t.log()-1/(12*t**2)
    assert minus>0
    delta_f = 6*t*(q(2).log()-3*p*minus)+upper
    assert delta_f < 0
    theta = osc["angle"](11431)
    radial = (q(1)/6).union(q(1)/4)
    gain = osc["phase_energy"](4*p+1,theta,radial)-osc["phase_energy"](p,theta,radial)
    assert gain + q("5e-18")*delta_f < 0
    for n in (2,5,14):
        product = q(1)
        for j in range(1,n):
            product *= 2*(j*q.pi()/(2*n)).sin()
        assert product.overlaps(q(n).sqrt())
    print("middle temperature cost per oscillator unit",lower.str(10),upper.str(10),
          "92905-step free energy change upper at t=1.1", (gain+q("5e-18")*delta_f).upper().str(12),flush=True)


if __name__=="__main__":
    verify()
