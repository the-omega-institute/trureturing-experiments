"""在相同的一圈处方下，不加有限径向质量项的亚普朗克联合算例。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

radial = run_path(str(Path(__file__).with_name("ch39-radial-fermion-feedback.py")))
free = radial["free"]


def verify():
    ctx.prec = 165
    r = arb("4e-5")
    relative = arb("3e-18")*(arb(25000)/42)**4
    h = arb("5e-18")
    u = h*relative
    epsilon2 = arb("5e-10")
    assert arb("3e-18") < 4*arb.pi()*epsilon2**2
    scale = 1000*arb("26.81676714125598").exp()
    assert scale*25000 < free["absolute"]["gauge"]["shortest"]["model"]["fixed"]["geo"]["BOOK_PLANCK"]
    p = 23226
    theta = free["absolute"]["angle"](11431)
    alpha = 2*arb.pi()/3
    delta = alpha-p*theta
    norm = 4*r+4*r*r
    first = 2*r*(1+2*r)
    remainder = norm**2/(2*(1-norm))
    coff = 4/(8*arb.pi()**2)
    slope_mass = 4*coff*r*r
    first_error = 3*coff*remainder*first
    second_error = 3*coff*((-2*(1-2*r).log())*first**2+remainder*2*r*r)
    d = arb(4).union(arb("4.8"))
    def radial_force(s):
        return 3*relative*(2*d*(s-1)+4*s*(s*s-1))+slope_mass*s+arb(0,first_error.upper())
    center_low,center_high = arb(".9999820"),arb(".9999838")
    assert radial_force(center_low)<0 and radial_force(center_high)>0
    half_curvature = (3*relative*(2*d.lower()-4)+slope_mass-second_error)/2
    assert half_curvature > 0
    # 两角能差为积分，其导数的绝对值上下界来自质量平方的全部本征值。
    def angle_slopes(radius):
        base = coff*r**3*radius**2
        largest = max(radius.upper(), arb(1))
        return (base/(1+r*largest)**2,base/(1-r*largest)**2)
    center = center_low.union(center_high)
    _, max_slope = angle_slopes(center)
    one_minus_cos = 2*(delta/2).sin()**2
    trial = (max_slope*one_minus_cos
             + 3*relative*((center**2-1)**2+2*one_minus_cos*center**2)/p
             + u*3*arb.pi()*(p-1)/2)
    width = (trial.upper()/half_curvature).sqrt()
    box = center+arb(0,width.upper())
    lower,upper = angle_slopes(box)
    force_early = 2*lower*(theta/2).sin()*min((delta+theta/2).sin().lower(),alpha.sin().lower())
    early = force_early-u*3*arb.pi()/2
    middle_distance = min((theta-delta).lower(),(2*alpha-(p+21000)*theta).lower())
    assert (1-middle_distance.cos())/(p+21000)>one_minus_cos/p
    budget = 3*relative*abs(box**2-1).upper()**2/p
    middle = u*3*arb(2).log()-budget/p
    # 远段可以失去全部相位能与至多预算/(p+21000) 的每步径向能。
    late = (u*3*arb(2).log()-budget/(p+21000)
            -(upper*one_minus_cos+6*relative*box.upper()**2*one_minus_cos/p)/21000)
    assert min(early.lower(),middle.lower(),late.lower()) > 0
    phase_shift = 12*r*abs(box-1).upper()/(1-r*max(box.upper(),arb(1)))
    each_shift = arb(4)/3*phase_shift/(16*arb.pi()**2)
    radial["fixed_gauge_box"](r,each_shift)
    assert r*42/(epsilon2.sqrt()*25000*box.lower()) < arb(".00301")
    print("对齐半径范围",center_low.str(8),center_high.str(8),
          "竞争者额外半宽",width.str(10),"长段每步能隙下界",late.lower().str(10),
          "三支额外原始阈值界",each_shift.str(10))


if __name__ == "__main__":
    verify()
