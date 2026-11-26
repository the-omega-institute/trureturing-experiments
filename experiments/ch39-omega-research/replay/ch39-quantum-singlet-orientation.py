"""规范单态量子热迹中的计次与非对易取向观测。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

quantum = run_path(str(Path(__file__).with_name("ch39-quantum-messenger-lengths.py")))


def verify():
    quantum["verify"]()
    ctx.prec = 160
    R = arb("3e-18")*(arb(25000)/42)**4
    u = arb("5e-18")*R
    hot = u/2
    r, C = arb("4e-5"),4/(8*arb.pi()**2)
    tau = 6*r+9*r**2
    rem = 3*C*tau**3/(6*(1-tau))
    edge = (u*42/25000)**2/4
    mediator = edge/arb("5e-10")
    p = 23226
    theta = quantum["classical"]["thermal"]["shape"]["scalar"]["free"]["absolute"]["angle"](11431)
    delta = 2*arb.pi()/3-p*theta
    end_trial = R*12*(delta/2).sin()**2/p
    guard = arb(1)/4
    # 保留第一条链端能，取出标量质量主项后剩余半份路径罚项非负。
    barrier = (R*((p*theta).sin()-guard)**2/p
               - 2*rem - 4*(C*r**2)**2/(4*R)-end_trial)
    assert (p*theta).sin() > guard and barrier > arb("5e-12")
    band, weight = arb("1e-12"),arb(200000)
    guard_margin = barrier-band-8*mediator*weight**2
    assert guard_margin > band and 2*guard*weight == 100000

    # 双温迹仅用于固定的目标长度。中性振子熵的损失是一份有限和。
    oscillator = 3*u*(p-1)*arb(p+1).log()
    packet = arb("1.51e-22")
    double = 2*hot
    v8 = arb.pi()**4/24
    log_heat = (9*v8.log()+24*arb(2).log()+12*(double/(4*arb.pi()*mediator)).log()
                +24*(double/(4*arb.pi()*edge)).log())
    edge_volume = 4**8+24*(double/R)**4
    high = (band-oscillator-2*packet)/double-log_heat-6*edge_volume.log()
    assert high > arb("5e11")
    print("Orientation barrier",barrier.lower().str(10),"low-state exponent",100000,
          "high-energy trace exponent",high.lower().str(10))


if __name__ == "__main__":
    verify()
