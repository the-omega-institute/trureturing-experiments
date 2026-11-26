"""三态费米子对中介共同半径的回馈与可变长度的自洽检验。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

free = run_path(str(Path(__file__).with_name("ch39-free-winding-cost.py")))


def verify():
    ctx.prec = 170
    r = arb(1)/2000
    relative = arb("1e-18")*(arb(10000)/42)**4
    h = arb("1e-17")
    u = h*relative
    p, q = 23226,92905
    theta = free["absolute"]["angle"](11431)
    alpha = 2*arb.pi()/3
    delta = alpha-p*theta
    rival = abs(q*theta-4*alpha)
    # 四个弱分量，矩阵以 M 为单位，0 <= s <= 2。
    coefficient = 4/(8*arb.pi()**2)
    norm_E = 4*r+4*r*r
    norm_first = 2*r*(1+2*r)
    norm_second = 2*r*r
    nonlinear = norm_E**2/(2*(1-norm_E))
    log_bound = -2*(1-2*r).log()
    first_error = 3*coefficient*nonlinear*norm_first
    second_error = 3*coefficient*(log_bound*norm_first**2+nonlinear*norm_second)
    mass_term = coefficient*4*r*r
    needed = (h*arb.pi()*p*q/2).sqrt()
    assert needed < arb(".000185")
    # 此窄窗中的树级径向力至多如下；费米子力为正且显著更大。
    radius = (1-needed).sqrt().union((1+needed).sqrt())
    d = arb(4).union(arb("4.8"))
    c = 1+arb(1)/p
    e = 2*(1-delta.cos())/p
    tree_force = 3*relative*(2*d*(radius-1)+2*e*radius+4*c*radius*(radius**2-1))
    minimum_force = mass_term*radius.lower()-first_error-abs(tree_force).upper()
    assert minimum_force > arb("4.99e-8")
    assert free["vacuum"](rival)<free["vacuum"](delta)
    # 一次共同的对中校准消去 c=1,s=1 处的费米子径向力。
    stiffness = 3*relative*d.lower()-(second_error+first_error)/2
    assert stiffness>arb("3.8e-8")
    static = free["vacuum"](delta)-free["vacuum"](arb(0))+3*relative*e
    maximum = static+u*(3*arb.pi()/2)*(p-1)
    window = (maximum/stiffness).sqrt()
    assert window < arb("4.0e-7")
    lower_ratio = (1-window)**2*((1-r)/(1+r*(1+window)))**6
    assert lower_ratio > arb(".99")
    loss = maximum/99
    lower,upper = free["cost_window"]()
    early = upper-u*3*arb.pi()/2
    cut = min((theta-delta).lower(),(2*alpha-(p+21000)*theta).lower())
    middle = free["vacuum"](cut)-free["vacuum"](delta)
    late = (u*3*arb(2).log()-lower)*21000-3*relative*e
    assert min(early.lower(),middle.lower(),late.lower()) > 10*loss
    opposite = free["vacuum"]((p*theta+alpha).cos().acos())-free["vacuum"](delta)
    assert opposite>10*loss
    # 改变被选半径时，单位模奇异值改变至多 r*window。
    # 对任意动量的三条对数核比较，再保守计入两个减法点。
    change = 12*r*window/(1-r*(1+window))
    each_shift = (arb(4)/3)*change/(16*arb.pi()**2)
    assert each_shift < arb("2.1e-11")
    fixed_gauge_box(r, each_shift)
    print("未补偿时的径向力下界",minimum_force.str(10),
          "振子必要的平方模容差",needed.str(10),
          "统一校准后的半径误差",window.str(11),
          "每个原始逆规范耦合的额外阈值上界",each_shift.str(11))


def fixed_gauge_box(r, each_shift):
    # 复算原尺度与编号的排除，加入新的阈值误差。
    absolute = free["absolute"]
    gauge = absolute["gauge"]
    absolute["ground"]["flow_bounds"]()
    time = arb("26.81676714125598")
    tol = arb(2)**-21
    model = gauge["shortest"]["model"]
    target = model["target"](11431,time)
    ratios = gauge["shortest"]["ratios"]
    base = gauge["gauge_state"](gauge["word"],time,(18,6,3))
    mixing = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    residuals = [1/(x*v)-target+arb(0,m.upper()) for x,v,m in zip(base,ratios,mixing)]
    rates = [arb(lo).union(arb(hi)) for lo,hi in (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    weights = [5*ratios[0],3*ratios[1],-8*ratios[2]]
    proj = sum(abs(w)*(abs(x).upper()+tol) for w,x in zip(weights,residuals))
    old_split = ((1+r**3)/(1-r**3)).log()/arb.pi()**2
    scale_width = (proj+old_split+16*each_shift)/abs(sum(w*z for w,z in zip(weights,rates))).lower()
    angles = arb(0).union(arb.pi())
    span = arb(0,scale_width.upper())
    shifts = absolute["shifts"](arb(42),span,time,angles,size=3,hop=r)
    color = target+residuals[2]+rates[2]*span-shifts[2]+arb(0,each_shift.upper())
    assert color+tol < model["target"](11430,time+scale_width)
    assert color-tol > model["target"](11432,time-scale_width)
    offset = arb("4.3e-7")
    shifts = absolute["shifts"](arb(42),offset,time,angles,size=3,hop=r)
    witness = [x+v*offset+target*(1-(-2*offset).exp())-s/g
               for x,v,s,g in zip(residuals,rates,shifts,ratios)]
    assert max(abs(x).upper() for x in witness)+each_shift/min(ratios)<tol
    # 固定此尺度时，两类三态在整个径向盒内均可通过原规范窗。
    mass_square = ((1-2*r)**2).union((1+2*r)**2)
    def summed_kernel(q_over_M):
        z = q_over_M**2/mass_square
        return 3*(z/5-3*z**2/140).union(z/5)
    common_box = (summed_kernel(offset.exp()/42)-summed_kernel((-time).exp()/42))/(12*arb.pi()**2)
    uniform_witness = [x+rate*offset+target*(1-(-2*offset).exp())-common_box/v
                       for x,rate,v in zip(residuals,rates,ratios)]
    assert all(abs(x)<tol for x in uniform_witness)



if __name__ == "__main__":
    verify()
