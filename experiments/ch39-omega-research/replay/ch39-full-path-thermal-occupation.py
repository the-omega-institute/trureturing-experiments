"""把三条路径的六个原始矩阵边也纳入无界经典配置积分。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

thermal = run_path(str(Path(__file__).with_name("ch39-matrix-endpoint-thermal-occupation.py")))


def verify():
    thermal["verify"]()
    ctx.prec = 160
    r, C = arb("4e-5"),4/(8*arb.pi()**2)
    R = arb("3e-18")*(arb(25000)/42)**4
    u, gap = arb("5e-18")*R,arb("3.5e-20")
    hot = u/2
    p, tail = 23226,10**16
    step = 3*u*arb(2).log()
    omega = thermal["shape"]["scalar"]["free"]["absolute"]["angle"](11431)
    delta = 2*arb.pi()/3-p*omega
    end = R*12*(delta/2).sin()**2/p
    tau = 6*r+9*r**2
    remainder = 3*C*tau**3/(6*(1-tau))
    # 外边模长大于四时，可留出 R y² 做高斯积分，余量仍超过试探。
    unit_cost = 2*C*r*r+2*remainder+end
    assert R*(4*(4**2-1)**2-4**2) > gap+unit_cost+(p-1)*3*u*arb.pi()/2
    assert step*(tail-p)-unit_cost > gap

    # 分母在最佳矩阵与达到原路径配方等号的六边旁积分。
    # 以下界只用原边及中介范数均不超过 2，余下的项由端点证书控制。
    sqrt2 = arb(2).sqrt()
    old_edges = 6*6*2*(1+4*sqrt2)*4*sqrt2
    attachments = 3*8*2*(2+4*sqrt2)*(1+4*sqrt2)
    assert old_edges+attachments+600 < 6000
    ball = arb("1e-18")
    edge_volume = 4**8 + 24*(hot/R)**4
    length_sum = tail + 1/(1-(-step/hot).exp())
    exponent = ((gap-6000*R*ball)/hot - 24*(2/ball).log() - 48*(1/ball).log()
                - 6*edge_volume.log() - length_sum.log())
    assert exponent > 31000
    print("Nine matrix variables, 72 real directions, wrong-sector lower exponent",
          exponent.lower().str(12))


if __name__ == "__main__":
    verify()
