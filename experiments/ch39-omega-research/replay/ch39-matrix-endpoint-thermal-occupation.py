"""对三个经典矩阵中介直接作热积分后，定量控制错误次数及反取向。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

shape = run_path(str(Path(__file__).with_name("ch39-general-messenger-alignment.py")))


def verify():
    gap = "3.5e-20"
    shape["verify"](gap)
    ctx.prec = 160
    gap = arb(gap)
    r, C = arb("4e-5"), 4/(8*arb.pi()**2)
    R = arb("3e-18")*(arb(25000)/42)**4
    u = arb("5e-18")*R
    tau = 6*r+9*r**2
    p = 23226
    theta = shape["scalar"]["free"]["absolute"]["angle"](11431)
    delta = 2*arb.pi()/3-p*theta
    end = R*12*(delta/2).sin()**2/p
    tail_start = 10**16
    step = 3*u*arb(2).log()
    rem = 3*C*tau**3/(6*(1-tau))
    assert step*(tail_start-p)-2*C*r**2-2*rem-end > gap
    # 每条矩阵以同一归一 Frobenius 半径移动时的总势变动。
    tree = 3*(4*arb("4.8")+24+64+40)+72
    loop = 8*C*r**2+3*C*tau**2/(2*(1-tau))*2*arb(2).sqrt()*r*(1+3*r)
    assert tree+loop/R < 600
    ball = arb("1e-17")
    hot = u/2
    prefactor = tail_start + 1/(1-(-step/hot).exp())
    exponent = (gap-600*R*ball)/hot - 24*(2/ball).log() - prefactor.log()
    assert exponent > 33000
    aligned_error = (gap/(arb(".48")*R)).sqrt()
    assert aligned_error < arb("4.41e-7")
    # 在能隙以内，原矩阵首角与中心投影首角的虚部同号。
    assert arb(".997")**2 * min(arb(2).sin(),arb("2.2").sin()) > aligned_error**2
    print("Lower exponent against wrong count or orientation",exponent.lower().str(11),
          "traceless norm inside the near-minimum window",aligned_error.upper().str(9))


if __name__ == "__main__":
    verify()
