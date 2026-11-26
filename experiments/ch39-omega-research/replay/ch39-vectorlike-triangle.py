"""三份矢量型场在同一有限参数点的壳相位、热权与绝对规范读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
absolute = run_path(str(root/"ch39-fermion-absolute-matching.py"))
gauge = absolute["gauge"]


def verify():
    ctx.prec = 180
    time = arb("26.81676714125598")
    scale = 1000*time.exp()
    planck = gauge["shortest"]["model"]["fixed"]["geo"]["BOOK_PLANCK"]
    upper_vev = (arb("1.001e-9")*4*arb.pi()).sqrt()
    five_field_minimum = 33*scale/upper_vev
    assert five_field_minimum/planck > arb("10.675")
    mass_ratio = 42
    hop = arb(1)/2000
    yukawa = hop*mass_ratio/((arb("5e-10").sqrt())*1000)
    assert yukawa/arb(".999") < 1 and 1000*scale < planck
    # 取中心扇区的无量纲中介模长，沿用其显式径向下界。
    power = 23226
    angle = absolute["angle"]
    eps = abs(power*angle(11431)-2*arb.pi()/3)
    next_eps = abs(power*angle(1144)-14*arb.pi()/3)
    assert 1/(1+(1-eps.cos())/(2*power)) > arb(".999")
    def energy(phase):
        squares = [1+hop**2+2*hop*((phase+2*k*arb.pi())/3).cos() for k in range(3)]
        return -4*sum(x**2*x.log() for x in squares)/(16*arb.pi()**2)
    cosine_gap = eps.cos()-next_eps.cos()
    contact_crossing = (energy(next_eps)-energy(eps))/cosine_gap
    determinant_lower = 1-2*hop**3
    determinant_upper = 1 + 3*hop**2/4 + 3*hop**4 + hop**6 + 2*hop**3
    derivative = 4*hop**3/(8*arb.pi()**2)
    assert derivative/determinant_upper < contact_crossing < derivative/determinant_lower
    assert arb("6.33257e-12") < contact_crossing < arb("6.33258e-12")
    cubic = contact_crossing*yukawa**2/(2*hop**3)
    assert arb(".02234") < cubic < cubic/arb(".999")**2 < arb(".02239")
    margin = arb("2e-6")
    kappa = arb("1e23")
    inverse_temperature = arb("2e22")
    tail = 22*arb.pi()*kappa.sqrt()*(arb(1)/2-2*kappa*(margin/(2*(power-1))).sin()**2).exp()
    weight_ratio = 2*(1+arb.const_euler()/2-(4*arb.pi()).log()/2)*(arb(1)/4+acb.zeta_zero(11431).imag**2)
    error = (weight_ratio/(1-tail)*
             ((-inverse_temperature*(energy(next_eps-margin)-energy(eps+margin))).exp()
              +tail*(inverse_temperature*(energy(eps+margin)-energy(arb(0)))).exp()))
    assert error < arb("4.12e-68")

    # 每支在同一质量谱与尺度上作欧氏减法，计入好角窗内全部热偏移。
    offset = arb("4.3e-7")
    phase = arb(0).union(eps+margin)
    shifts = absolute["shifts"](arb(mass_ratio),offset,time,phase,size=3,hop=hop)
    base = gauge["gauge_state"](gauge["word"],time,(18,6,3))
    target = gauge["shortest"]["model"]["target"](11431,time)
    ratios = gauge["shortest"]["ratios"]
    mixing = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    absolute["ground"]["flow_bounds"]()
    speeds = [arb(lo).union(arb(hi)) for lo,hi in
              (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    boxes = [1/(v*x)-target*(-2*offset).exp()+speed*offset+arb(0,m.upper())-shift/v
             for v,x,speed,m,shift in zip(ratios,base,speeds,mixing,shifts)]
    worst = max(abs(x).upper() for x in boxes)
    assert worst < arb(2)**-21
    print("五态质量所需原边标尺 / 书中普朗克标尺下界",(five_field_minimum/planck).str(11),
          "三态算例的带荷 Yukawa 上界",(yukawa/arb(".999")).str(10),"联合失败率",
          error.str(11),"好角窗内三支绝对残差上界",worst.str(11),
          "反转两个中心能级的无量纲接触系数",contact_crossing.str(12))


if __name__ == "__main__":
    verify()
