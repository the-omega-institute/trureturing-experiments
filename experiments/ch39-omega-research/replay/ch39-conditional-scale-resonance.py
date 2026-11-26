"""用标量条件热更新从全尺度区间准备壳上的首个共振。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

source = run_path(str(Path(__file__).with_name("ch39-moving-scale-resonance.py")))
q = arb


def verify():
    ctx.prec = 112
    joint = source["joint"]
    first, selected, last = (source[k] for k in ("first", "selected", "last"))
    base = selected-first
    model = joint["flow"]["gauge"]["shortest"]["model"]
    time, eps = q("26.81676714125598"), q("1.005e-7")
    reference = model["target"](selected,time)
    offset = model["fixed"]["OFFSET"]
    offset = q(offset.numerator)/offset.denominator
    heat_min, heat_max = q("1.25e-18"), q("2e-17")
    heights = [z.imag for z in acb.zeta_zeros(first,last-first+1)]
    assert len(heights) == 1278
    denominators = [h*h+offset for h in heights]
    shifts = [(denominators[base]/v).log()/2 for v in denominators]
    angles = [2*(1/(2*h)).atan() for h in heights]
    alpha = 2*q.pi()/3
    counts = [int((alpha/a+q(1)/2).floor().unique_fmpz()) for a in angles]
    phase_ceiling = [6*(1-(p*a-alpha).cos())/p for a,p in zip(angles,counts,strict=True)]
    speeds = [q(lo).union(q(hi)) for lo,hi in
              (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    mean = sum(speeds)/3
    pairs = [speeds[j]-speeds[i] for i,j in combinations(range(3),2)]
    assert mean > 0 and all(v > 0 for v in pairs)
    largest_shift = max(abs(b).upper() for b in shifts)
    small = q("1e-10")
    band = q(".01")
    assert all(23 < time+b-band and time+b+band < 33 for b in shifts)
    residual = lambda w,b: eps+mean.upper()*(b+w)+reference*((2*w).exp()-1)
    contrast = lambda w,b: sum((2*eps+v.upper()*(b+w))**2 for v in pairs)/3
    exterior = lambda w,b: reference*(1-(-2*w).exp())-eps-mean.upper()*(b+w)
    # 实际定位只需各编号自己的 Q 和三支规范数据。先用单调残差求零点。
    assert exterior(q(".001"),largest_shift) > 0
    assert q(10)/2**24 < q(".001")
    for count,k in enumerate((".023",".026",".030"),1):
        m,d = (q(v) for v in joint["leads"]["terminal_modes"](count))
        widths = [q("2e-7")+q(k)*abs(b) for b in shifts]
        tails, margins = [], []
        for b,w in zip(shifts,widths,strict=True):
            assert small < w < band/2
            near = 3*m*residual(small,abs(b))**2 + d*contrast(small,abs(b))
            r = exterior(w,abs(b))
            assert r > 0
            tails.append(((3*m*r*r-near)/heat_max-(5/small).log()).lower())
        for j in range(len(heights)-1):
            inside,outside = (j+1,j) if j<base else (j,j+1)
            A,B = abs(shifts[inside]),abs(shifts[outside])
            w = widths[outside]
            def square_drop(v, error):
                return (B-A)*(v.lower()**2*(B+A)-2*v.upper()*error)
            common_error = eps + mean.upper()*w + reference*((2*w).exp()-1)
            drop = 3*m*square_drop(mean,common_error)
            drop += d*sum(square_drop(v,2*eps+v.upper()*w) for v in pairs)/3
            prior = heat_max*((heights[j+1]**2+q(1)/4)/(heights[j]**2+q(1)/4)).log()
            margin = drop-phase_ceiling[inside]-abs(counts[j+1]-counts[j])*q("2e-16")-prior
            margins.append((margin/heat_max).lower())
        assert min(tails)>23000 and min(margins)>20000
        print("引线",count,"尺度带外负对数下界",min(tails).str(9),
              "共振边向内能降下界",min(margins).str(9))

        # 在半宽 0.01 内严格凸。带外的共同端口项在两侧单调远离零。
        r_low = exterior(band/2,largest_shift)
        center_top = 3*m*(eps+mean.upper()*largest_shift)**2+d*contrast(0,largest_shift)
        assert r_low > 0 and 3*m*r_low*r_low-center_top > q("1e-4")
        q_span = reference*(-2*band).exp().union((2*band).exp())
        slope = 2*q_span+mean
        acc = 1+4*q_span.upper()
        curvature = 6*m*(slope**2+q(0,(residual(band,largest_shift)*acc).upper()))
        for v in pairs:
            pair_size = 2*eps+v.upper()*(largest_shift+band)
            curvature += 2*d*(v**2+q(0,(2*pair_size).upper()))/3
        assert curvature>20 and curvature<50

    # 拒绝采样：99% 在近极小点的半宽 2 sqrt(kBT) 带内，1% 为全域均匀。
    center_error = q("1e-11")
    assert band/2**32 < center_error
    envelope_per_sqrt_heat = 4/q(".99")
    mass_per_sqrt_heat = 2/q(50).sqrt()*(-q(1)/2-25*center_error**2/heat_min).exp()
    assert envelope_per_sqrt_heat / mass_per_sqrt_heat < 24
    outside_weight = (-10*(2-center_error/heat_min.sqrt())**2).exp()
    assert q("1e-4")/heat_max > 40
    assert 1000*outside_weight < envelope_per_sqrt_heat*heat_min.sqrt()
    # 每次拒绝采样成功率大于 1/24，最多五百万次试探的总时限。
    cycles = 180000
    expected = q(5000000)/24
    budget_exponent = (expected-cycles)**2/(2*expected)
    assert budget_exponent > 1900
    growth = 2*(q(160000)/4-max(counts)+1)**2/160000
    drift = 2*(q(20000)/4-max(selected-first,last-selected))**2/20000
    failure = ((-growth).exp()+(-drift).exp()+(-budget_exponent).exp()
               +cycles*((-q(23000)).exp()+(-q(20000)).exp()+(-q(64)).exp()))
    assert failure < q("3e-23")
    print("任意初始尺度的失败概率上界",failure.upper().str(10))
    return shifts, reference, mean, pairs


if __name__ == "__main__":
    verify()
