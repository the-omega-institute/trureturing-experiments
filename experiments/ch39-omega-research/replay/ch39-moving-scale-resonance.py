"""连续规范尺度的错号陷阱及沿首个相位共振的联合制备。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

previous = run_path(str(Path(__file__).with_name("ch39-interleaved-winding-updates.py")))
joint = previous["joint"]
first, selected, last = previous["first"], previous["selected"], previous["last"]
q = arb


def verify():
    ctx.prec = 128
    model = joint["flow"]["gauge"]["shortest"]["model"]
    nominal = q("26.81676714125598")
    eps, max_heat = q("1.005e-7"), q("2e-17")
    target = model["target"]
    reference = target(selected, nominal)
    leads = [tuple(q(v) for v in joint["leads"]["terminal_modes"](k)) for k in (1, 2, 3)]
    offset = model["fixed"]["OFFSET"]
    offset = q(offset.numerator) / offset.denominator
    heights = [z.imag for z in acb.zeta_zeros(first, last - first + 1)]
    assert len(heights) == last - first + 1
    base = selected - first
    denominator = [h*h + offset for h in heights]
    angles = [2*(1/(2*h)).atan() for h in heights]
    shifts = [(denominator[base]/d).log()/2 for d in denominator]
    alpha = 2*q.pi()/3
    counts = [int((alpha/a + q(1)/2).floor().unique_fmpz()) for a in angles]
    assert counts[base] == 23226 and (min(counts), max(counts)) == (22139, 24386)
    assert min(b-a for a,b in zip(counts, counts[1:])) >= 0
    four_step_edges = [first+j for j in range(len(counts)-1) if counts[j+1]-counts[j] == 4]
    assert max(b-a for a,b in zip(counts, counts[1:])) == 4 and four_step_edges

    # 固定长度时给 11430 一个保留算术规范量的试探尺度。
    angle_diff = angles[base-1] - angles[base]
    min_angle = alpha - 23000*angles[base-1]
    assert 0 < min_angle < alpha
    phase_saving = (6*angle_diff*min(alpha.sin().lower(), min_angle.sin().lower())
                    / (1+q(3)/8)**2)
    trial_gauge = 3*max(d for _,d in leads)*(eps+q(".09")*abs(shifts[base-1]))**2
    assert (phase_saving-trial_gauge)/max_heat > 16000000

    # 11432 初态小带相对其外部的尺度权。编号提议若伴随平移，仍计入其最坏误差。
    trap = base + 1
    shift = shifts[trap]
    inner, outer, move = q("1e-8"), q("6e-6"), q("1e-6")
    local_error = lambda w: eps+q(".09")*(abs(shift)+w)+reference*((2*w).exp()-1)
    tail = reference*(1-(-2*outer).exp())-eps-q(".09")*(abs(shift)+outer)
    assert tail > 0
    for m,d in leads:
        exponent = (3*m*tail**2-3*d*local_error(inner)**2)/max_heat-(5/inner).log()
        assert exponent > 18000000
        for other in (trap-1,trap+1):
            other_q = reference*denominator[trap]/denominator[other]
            step = abs(other_q-reference)*(-2*outer).exp()
            moved = other_q*(2*outer).exp()*((2*move).exp()-1)+q(".09")*move
            after = step-local_error(outer)-moved
            assert after > 0
            contrast = 3*(eps+q(".09")*(abs(shift)+outer))**2
            gauge_gain = 3*m*(after**2-local_error(outer)**2)-d*contrast
            prior = max_heat*abs(((heights[trap]**2+q(1)/4)/(heights[other]**2+q(1)/4)).log())
            assert (gauge_gain-6*abs(angles[trap]-angles[other])-prior)/max_heat > 100000000

    # 以已有的三支规范速度带估计保持 Q 的联合动作。
    bands = [q(a).union(q(b)) for a,b in
             (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    mean_speed = sum(bands)/3
    differences = [bands[j]-bands[i] for i,j in combinations(range(3),2)]
    assert mean_speed > 0 and all(v > 0 for v in differences)
    measurement = q("5e-7")
    gauge_radius = 2*eps/q(".1435")
    tolerance = eps+q(".09")*gauge_radius+measurement
    zmax = -(1-tolerance/reference).log()/2
    # 各编号移步保持 Q。独立尺度步即使连续同向，对有限预算的最坏影响仍在校准窗内。
    initial_accuracy = q("2e-7")
    independent_step = q("2e-13")
    total_steps = 180000
    drift_in_q = (reference+tolerance)*((2*total_steps*independent_step).exp()-1)
    assert initial_accuracy+drift_in_q < measurement
    assert all(23 < nominal+s-zmax and nominal+s+zmax < 33 for s in shifts)
    mean_error = eps+q(".09")*zmax+tolerance
    pair_error = 2*(eps+q(".09")*zmax)

    # 逐一控制首个共振前的缩短与共振处的离开。
    phase_tops = []
    length_margins = []
    local_margins = []
    radial = (q(1)/6).union(q(1)/4)
    for a,p in zip(angles, counts, strict=True):
        delta = alpha-p*a
        assert abs(delta) < a/2
        u = lambda k: 2*(1-(alpha-k*a).cos())/k
        middle, lower, upper = u(p), u(p-1), u(p+1)
        phase_tops.append((3*middle).upper())
        energy_change = lambda v: 3*(v-middle)/((1+radial*v/2)*(1+radial*middle/2))
        # 除末步外，最靠近共振的中点角距为 delta+3a/2。
        earlier = (12*(a/2).sin()*min(alpha.sin().lower(), (delta+3*a/2).sin().lower())
                   / (p*(1+q(3)/8)**2))
        length_margins.extend([((earlier-q("2e-16"))/max_heat).lower(),
                               ((energy_change(lower)-q("2e-16"))/max_heat).lower(),
                               ((energy_change(upper)+q("1e-17"))/max_heat).lower()])
        # 从一端半径四格找另一端，最远候选可高出另一端首个共振八格。
        assert 0 < (p+8)*a-alpha < q.pi()/3
        old = upper
        for k in range(p+2,p+9):
            new = u(k)
            rise = 3*(new-old)/((1+radial*old/2)*(1+radial*new/2))
            local_margins.append(((rise+q("1e-17"))/max_heat).lower())
            old = new
    assert min(length_margins) > 64
    assert min(local_margins) > 64

    label_margins = []
    for j in range(len(heights)-1):
        inside, outside = (j+1,j) if j < base else (j,j+1)
        small, large = abs(shifts[inside]), abs(shifts[outside])
        step = large-small
        assert step > 0
        def square_drop(speed, error):
            return step*(speed.lower()**2*(large+small)-2*error*speed.upper())
        prior_change = max_heat*((heights[j+1]**2+q(1)/4)/(heights[j]**2+q(1)/4)).log()
        deductions = phase_tops[inside]+abs(counts[j+1]-counts[j])*q("2e-16")+prior_change
        for m,d in leads:
            gauge_drop = (3*m*square_drop(mean_speed,mean_error)
                          +d*sum(square_drop(v,pair_error) for v in differences)/3)
            label_margins.append(((gauge_drop-deductions)/max_heat).lower())
    assert min(label_margins) > 100000

    growth_steps, label_steps = 160000, 20000
    assert growth_steps + label_steps == total_steps
    growth = 2*(q(growth_steps)/5-max(counts)+1)**2/growth_steps
    drift = 2*(q(label_steps)/5-max(selected-first,last-selected))**2/label_steps
    risk = ((-growth).exp()+(-drift).exp()
            +(growth_steps+label_steps)*((-q(64)).exp()+(-q(100000)).exp()))
    assert risk < q("3e-23")
    print("11430 在 23000 以内比所有 11431 尺度至少低（热能单位）",
          ((phase_saving-trial_gauge)/max_heat).lower().str(10))
    print("长度边与联合编号边的最小能差下界（热能单位）",
          min(length_margins).str(10),min(label_margins).str(10))
    print("180000 步的失败概率上界",risk.upper().str(10))
    print("需要四格连接的边数",len(four_step_edges),"首条边的左端编号",four_step_edges[0],
          "共振后第二至第八步的最小涨幅（热能单位）",min(local_margins).str(10))


if __name__ == "__main__":
    verify()
