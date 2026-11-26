"""十四个弱顶点的 SU(2) 单态热迹，使用两奇异值径向势。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
q = arb


def verify(*, singlet_clocks=False):
    ctx.prec = 128
    heat_low, heat = q("1.25e-18"),q("5e-18")
    h_low, h_high = q("5e-18"),q("1e-17")
    k_max = q(10)
    length = 23226
    previous = run_path(str(root/"ch39-thermal-plaquette-occupation.py"))
    joint = previous["joint"]
    theta = joint["winding"]["angle"](11431)
    alpha = 2*q.pi()/3
    delta = alpha-length*theta
    chord = lambda x: 2*(1-x.cos())
    e = chord(delta)/length
    reverse = chord(alpha+length*theta)/length
    assert 0 < delta < theta/2
    components = 1 if singlet_clocks else 4
    assert components*(3*q.pi()/2)*h_high < q("2e-16")
    assert components*q("2.0794")*h_low > q("1e-17")
    slope = 3/(1+q(3)/8)**2
    early = slope*4*(theta/2).sin()*min((delta+theta/2).sin().lower(),alpha.sin().lower())/length-q("2e-16")
    middle = slope*(chord(theta-delta)/(length+21000)-e)
    tail = q(21000)*q("1e-17")-3*e
    assert theta-delta < 2*alpha-(length+21000)*theta
    count_power = min(v.lower() for v in (early,middle,tail))/heat-q(21010).log()
    assert early/heat>1 and 1/(1-(-q("1e-17")/heat).exp())<3

    time, width = q("26.81676714125598"),q("1e-8")
    separation_radius, eps = q(".00005"),q("1.005e-7")
    model = joint["flow"]["gauge"]["shortest"]["model"]
    target = model["target"]
    center = target(11431,time)
    local_error = eps+(q(".09")+2*center*(2*width).exp())*width
    separation = min((target(11430,time+separation_radius)-center).lower(),
                     (center-target(11432,time-separation_radius)).lower())-eps-q(".09")*separation_radius
    split = q(".1435")*separation_radius-2*eps
    assert separation>0 and split>0
    label_powers = []
    for j in (1,2,3):
        m,d = (q(v) for v in joint["leads"]["terminal_modes"](j))
        gauge = min((3*m*separation**2).lower(),(d*split**2/2).lower())
        trial = 3*e+length*q("2e-16")+3*d*local_error**2
        label_powers.append(((gauge-trial)/heat-(15/width).log()-16).lower())

    # 60 个复坐标的中心小盒，及归一矩阵高斯包给出的量子下界。
    box = q(".1")
    coord = 2*box
    field_sq = (1+coord*heat.sqrt())**2
    assert 2*coord+q(2).sqrt()*coord**2*heat.sqrt() < q(".41")
    box_energy = q(".41")**2*(6*q("7.5")+9*6+9*(k_max+q(4)/3)+12)
    path_cost = k_max + (30 if singlet_clocks else 2)
    linear = 15*6*(16*field_sq-8)+path_cost*24*field_sq
    quadratic = 15*6*32+path_cost*48
    packet_energy = 4*(150+linear/10+quadratic*h_high/100)
    span = 15*q(2).log()-60*(box**2*heat_low).log()+box_energy+packet_energy

    # 对所有十四个顶点取 Haar 平均。用平均能量估计每个单位李代数方向的绝对荷。
    energy = 3*e+2*heat*(span+60*q(2).log()+6*components*(length-1)*(2*q(2).cosh()).log())
    radius_sq = q(45)/2+energy/8
    if not singlet_clocks:
        radius_sq += 12*length+(q(length)**2+q(length)/2)*energy
    charge = (4*energy*radius_sq/h_low**2).sqrt()
    angle = 1/(28*charge)
    assert 0 < angle < q.pi()/2
    ball = 8*angle**3/(3*q.pi()**3)
    loss = q(2).log()-14*ball.log()
    count_gap = count_power-span-loss
    label_gap = min(label_powers)-span-loss
    chi_slope = 3*(reverse-e)/(2*(1+e/8)*(1+reverse/8))
    chi_error = 2*heat*(span+loss)/chi_slope
    print("量子体积对数跨度",span.str(10),"十四点投影损失上界",loss.str(10),
          "错计次、错编号的指数下界",count_gap.str(11),label_gap.str(11),
          "取向缺额",chi_error.str(10))
    assert count_gap>14000 and label_gap>790000 and chi_error<q("3e-10")
    return count_power,min(label_powers),e,reverse


if __name__ == "__main__":
    verify()
