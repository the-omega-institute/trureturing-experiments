"""在三个复标量中心背景上同时放开模长和相对相位。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

radial = run_path(str(Path(__file__).with_name("ch39-radial-fermion-feedback.py")))
free = radial["free"]


def verify(energy_slack="0", *, label=11431, count=23226, gauge_check=None):
    ctx.prec = 160
    slack = arb(energy_slack)
    assert slack >= 0
    r = arb("4e-5")
    energy = arb("3e-18") * (arb(25000)/42)**4
    u = arb("5e-18") * energy
    closure = 2
    coefficient = 4/(8*arb.pi()**2)
    m = 2*coefficient*r*r
    perturbation = 4*r+4*r*r
    remainder = perturbation**2/(2*(1-perturbation))
    first = 2*r*(1+2*r)
    first_error = 3*coefficient*remainder*first
    hess_error = 3*coefficient*((-2*(1-2*r).log())*first**2 + remainder*2*r**2)
    p = count
    theta = free["absolute"]["angle"](label)
    hot, cold = 3*arb(2).log()*u, 3*arb.pi()*u/2
    chord = lambda x: 4*(x/2).sin()**2
    margins = []
    all_angles = []
    all_radii = []
    for cell in range(16):
        d = (4 + arb(cell)/20).union(4 + arb(cell+1)/20)
        approx_shift = m/(energy*(2*d+8))
        residual = (first_error + arb(2).sqrt()*(m*approx_shift
                    + energy*(12*approx_shift**2+4*approx_shift**3)))
        weak = energy*(2*d.lower()-4)-hess_error
        strong = energy*(2*d.lower()+12*arb(".9999")**2-4)-hess_error
        assert 1-approx_shift-residual/weak > arb(".9999")
        a = 1-approx_shift+arb(0,(residual/strong).upper())
        # U 的极小点为 (a,a,1)。精确闭合的三角形用正取向试探。
        beta = (1/(2*a*a)-1).acos()
        gamma = arb.pi()-beta/2
        products = [a*a, a, a]
        directions = [beta, gamma, gamma]
        high = coefficient*r**3*a.upper()**2/(1-r)**2
        trial_static = high*chord(beta-p*theta)/2 + energy*sum(
            abs(h-1).upper()**2 + h.upper()*chord(b-p*theta)
            for h,b in zip(products,directions))/p
        budget = trial_static.upper() + cold*(p-1) + slack
        assert trial_static < arb('1.4e-19') and budget < arb('4e-19')
        assert a.lower()-(2*budget/weak).sqrt() > arb(".9999")
        width = (2*budget/strong).sqrt()
        z = [a+arb(0,width.upper()),a+arb(0,width.upper()),arb(1, width.upper())]
        assert all(s > arb('.9999') and s < arb('1.0001') for s in z)
        all_radii.append(z)
        gap = (budget/(closure*energy)).sqrt()
        # 对每条边先用闭合残差控制对边的长度，再用余弦定理。
        angles, lengths = [], []
        for i in range(3):
            x,y,w = z[i],z[(i+1)%3],z[(i+2)%3]
            cosine = ((w+arb(0,gap.upper()))**2-x*x-y*y)/(2*x*y)
            angles.append(cosine.acos())
            lengths.append(x*y)
        assert all(2<angle and angle<arb("2.2") for angle in angles)
        all_angles.append(angles[0])
        # 若能竞争，每个相位都离 q theta 的某个 2pi 平移不足 .001。
        max_q = 69000 if slack > 0 else 60000
        assert hot*(max_q-p)>trial_static + slack
        assert energy*min(h.lower() for h in lengths)*chord(arb(".001")) > max_q*budget
        windows = []
        for orientation in (1,-1):
            for turn in range(-1,3):
                location = (orientation*angles[0]+2*arb.pi()*turn+arb(0,arb(".001")))/theta
                lo = max(1, int(location.lower().floor().unique_fmpz()))
                hi = min(max_q-1, int(location.upper().ceil().unique_fmpz()))
                windows.extend((q,orientation) for q in range(lo,hi+1))
        assert max_q*theta < 2*arb.pi() and (p,1) in windows
        low = coefficient*r**3*lengths[0].lower()/(1+r*max(s.upper() for s in z))**2
        radial_budget = sum(max(arb(0),abs(h-1).lower())**2 for h in lengths)
        for q,orientation in windows:
            if (q,orientation)==(p,1):
                continue
            angle_budget = sum(h.lower()*max(arb(0),chord(orientation*b-q*theta).lower())
                               for h,b in zip(lengths,angles))
            fermion = low*max(arb(0),chord(orientation*angles[0]-q*theta).lower())/2
            osc = (hot if q>=p else cold)*(q-p)
            gap_q = fermion+energy*(radial_budget+angle_budget)/q+osc-trial_static.upper()
            assert gap_q>slack, (cell,q,orientation,gap_q)
            margins.append((gap_q.lower(),q,orientation))
    deviation = max(abs(z[i]-1).upper() for z in all_radii for i in (0,1))
    threshold = (arb(4)/3)*12*r*deviation/(1-r*(1+deviation))/(16*arb.pi()**2)
    if gauge_check is None:
        assert label == 11431
        radial["fixed_gauge_box"](r,threshold)
    else:
        gauge_check(r,threshold)
    angle_low = min(x.lower() for x in all_angles)
    angle_high = max(x.upper() for x in all_angles)
    print("Positive charged angle offset from 2pi/3", (angle_low-2*arb.pi()/3).str(9),
          (angle_high-2*arb.pi()/3).str(9))
    print("Smallest enumerated lower margin",min(margins)[0].str(9),"at",min(margins)[1:])


if __name__ == "__main__":
    verify()
