"""三个可独立变化的中心模长下，复算三态环的自由计次与规范筛选。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

radial = run_path(str(Path(__file__).with_name("ch39-radial-fermion-feedback.py")))
free = radial["free"]


def verify():
    ctx.prec = 165
    r = arb("4e-5")
    scale = arb("3e-18") * (arb(25000)/42)**4
    u = arb("5e-18") * scale
    coefficient = 4 / (8 * arb.pi()**2)
    mass_slope = 2 * coefficient * r**2
    d = arb(4).union(arb("4.8"))
    closure = 2
    a = 2*d + 8
    charged_shift = mass_slope/scale * (2/(3*a) + 1/(3*(a+3*closure)))
    third_shift = mass_slope/scale * (2*closure/(a*(a+3*closure)))
    # 全立方体、任意单位方向的费米子导数误差。
    perturbation = 4*r + 4*r**2
    remainder = perturbation**2 / (2*(1-perturbation))
    first = 2*r*(1+2*r)
    grad_error = 3*coefficient*remainder*first
    hess_error = 3*coefficient*((-2*(1-2*r).log())*first**2 + remainder*2*r**2)
    curvature = scale*(2*d.lower()-4) - hess_error
    assert curvature > 0
    v = charged_shift.upper()
    residual = (grad_error + arb(2).sqrt()*mass_slope*v
                + arb(3).sqrt()*scale*(12*v**2+4*v**3))
    anchor_error = residual/curvature
    centers = [1-z+arb(0,anchor_error.upper()) for z in (charged_shift,charged_shift,third_shift)]

    p = 23226
    theta = free["absolute"]["angle"](11431)
    alpha = 2*arb.pi()/3
    delta = alpha-p*theta
    chord = 2*(delta/2).sin()**2
    def budgets(box):
        prod = [box[i]*box[(i+1)%3] for i in range(3)]
        norm_hop = r*max(arb(1), *(x.upper() for x in box))
        slope = coefficient*r**3*prod[0]
        return (slope.lower()/(1+norm_hop)**2, slope.upper()/(1-norm_hop)**2,
                sum(abs(x-1).upper()**2 for x in prod), sum(x.upper() for x in prod))
    _, high, radial_cost, angle_cost = budgets(centers)
    trial = high*chord + scale*(radial_cost+2*chord*angle_cost)/p + u*3*arb.pi()*(p-1)/2
    width = (2*trial.upper()/curvature).sqrt()
    box = [s+arb(0,width.upper()) for s in centers]
    assert all(0<s and s<2 for s in box) and box[2] > box[0] and box[2] > box[1]
    low, high, radial_cost, angle_cost = budgets(box)
    low_angle = min((delta+theta/2).sin().lower(),alpha.sin().lower())
    early = 2*low*(theta/2).sin()*low_angle - u*3*arb.pi()/2
    cut = min((theta-delta).lower(),(2*alpha-(p+21000)*theta).lower())
    assert (1-cut.cos())/(p+21000) > chord/p
    middle = u*3*arb(2).log() - scale*radial_cost/p**2
    late = (u*3*arb(2).log() - scale*radial_cost/(p*(p+21000))
            -(high*chord+2*scale*angle_cost*chord/p)/21000)
    assert min(early.lower(), middle.lower(), late.lower()) > 0
    # 同一长度的两取向由端点的余弦角距严格排序。
    assert (p*theta+alpha).cos() < delta.cos()
    shift = max(abs(x-1).upper() for x in box[:2])
    each_shift = (arb(4)/3)*12*r*shift/(1-r*(1+shift))/(16*arb.pi()**2)
    radial["fixed_gauge_box"](r, each_shift)
    assert box[2].lower() - max(box[i].upper() for i in (0,1)) > arb("1.5e-5")
    print("charged radial ranges", box[0].lower().str(10), box[0].upper().str(10),
          "third range", box[2].lower().str(10), box[2].upper().str(10))
    print("anchor error", anchor_error.str(8), "competitive width", width.str(8),
          "positive per-step margins", early.lower().str(7),middle.lower().str(7),late.lower().str(7))


if __name__ == "__main__":
    verify()
