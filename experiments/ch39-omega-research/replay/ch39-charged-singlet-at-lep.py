"""关闭直接 Z 分宽的六种场表共用的带电标量在 LEP2 的截面与寿命尺度。"""
from flint import acb, arb, ctx


def production(energy, mass, alpha):
    z, width, sw2 = arb("91.1876"), arb("2.4952"), arb(".23129")
    s = energy**2
    beta = (1-4*mass**2/s).sqrt()
    # 取 S⁻ 为电荷 −1、T3=0，光子振幅归一成 1。
    propagator = s/acb(s-z**2, z*width)
    left = 1 + (sw2-arb(1)/2)/(1-sw2)*propagator
    right = 1 + sw2/(1-sw2)*propagator
    photon = arb.pi()*alpha**2*beta**3/(3*s)*arb("3.8937937e8")
    return photon, photon*(abs(left)**2+abs(right)**2)/2


def verify():
    ctx.prec = 192
    mass = arb("80.38")
    photon, total = production(arb("206.6"), mass, 1/arb(128))
    print("206.6 GeV 光子、光子与 Z 交换的 Born 截面/pb", photon.str(9), total.str(9))
    assert total > 5*arb(".028")
    # 即使使用零动量的电磁耦合，给定树级产生截面仍远大于 OPAL 上限。
    _, lower = production(arb("206.6"), mass, 1/arb("137.036"))
    assert lower > 4*arb(".028")
    hbar = arb("6.582119569e-25")  # GeV s
    for lifetime in (arb("1e-6"), arb("1e-3")/299792458):
        norm = (8*arb.pi()*hbar/(mass*lifetime)).sqrt()
        print("寿命/s 与质量忽略下的反对称双轻子耦合范数", lifetime.str(9), norm.str(9))
    # 每个反对称味对的两种末态中，至少一半的分宽给 e 或 μ。
    # 允许实测带电轻子质量，τ 的相空间不大于它的伙伴。
    weights = [(1-(arb(m)/mass)**2)**2 for m in (".000511",".10566","1.77686")]
    for i,wi in enumerate(weights):
        for j,wj in enumerate(weights[:i]):
            emu = (wi if i<2 else 0) + wj
            assert 2*emu > wi+wj
    print("纯双轻子衰变时同味 ee 或 μμ 的 Born 截面下界/pb", (total/8).str(8))


if __name__ == "__main__":
    verify()
