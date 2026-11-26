"""由正则端点动能和未消去的谐振链控制全正整数长度的量子热权。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

classical = run_path(str(Path(__file__).with_name("ch39-full-path-thermal-occupation.py")))


def verify():
    classical["verify"]()
    ctx.prec = 160
    R = arb("3e-18") * (arb(25000)/42)**4
    u = arb("5e-18")*R
    hot = u/2
    # 取 hbar Omega=M，原边与中介的正则动能标尺为 f_e、epsilon f_e。
    edge = (u*42/25000)**2/4
    mediator = edge/arb("5e-10")
    packet = arb("6e-11")
    root2 = arb(2).sqrt()
    # 下列值均为沿单位总平方范数端点变分的二阶 Taylor 系数界。
    radial = 6*(32+(1+4*root2)*2*root2)
    coupling = 8*(17+(2+4*root2)*root2)
    r = arb("4e-5")
    C = 4/(8*arb.pi()**2)
    tau = 6*r+9*r**2
    s1 = 2*root2*r*(1+3*r)
    s2 = 4*r*r
    loop = C*r*r + 3*C*((-2*(1-3*r).log())*s1**2
                         + tau**2/(2*(1-tau))*s2)/2
    assert radial+coupling+6+22+24+26+16+loop/R < 700
    potential = 700*R*72*packet**2
    kinetic = (24*mediator+48*edge)*arb.pi()**2/(4*packet**2)
    packet_cost = potential+kinetic
    gap = arb("3.5e-20")
    assert packet_cost < arb("1.6e-22")
    volume8 = arb.pi()**4/24
    log_heat_kernel = (9*volume8.log()+24*arb(2).log()
                       +12*(hot/(4*arb.pi()*mediator)).log()
                       +24*(hot/(4*arb.pi()*edge)).log())
    tail = 10**16 + 1/(1-(-3*u*arb(2).log()/hot).exp())
    edge_volume = 4**8 + 24*(hot/R)**4
    exponent = (gap-packet_cost)/hot - log_heat_kernel - 6*edge_volume.log()-tail.log()
    assert exponent > 34000
    print("Quantum length error exponent at least",exponent.lower().str(12),
          "packet energy overhead",packet_cost.upper().str(10))


if __name__ == "__main__":
    verify()
