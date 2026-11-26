"""三个物理窗口与一次中性壳响应共用五项算术时间的尺度下界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

neutral = run_path(str(Path(__file__).with_name("ch39-neutral-clock-response.py")))
model = neutral["model"]


def envelope(bound, left=0):
    ctx.prec, ctx.cap = 1024, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    moments = neutral["contrast"]["moments_five"](times)
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    left = arb(left)
    nodes, weights = neutral["fixed_ends"](moments, left, bound)
    spectrum = neutral["neutral_polynomial"]()
    spectrum = spectrum*(1/spectrum[0])
    # 在三个小节点匹配单位辐边的味导数，在大端点放宽到 19g。
    def score(v):
        return neutral["efficiency"](v, spectrum)*model["gradient"](v, poles, slopes)
    margin = sum(w*score(v) for w, v in zip(weights[:3], nodes[:3], strict=True))
    margin += 19*weights[3]*model["gradient"](bound, poles, slopes)
    assert margin < 0
    x = arb_poly([0, 1])
    factors = [1+v*x for v in poles]
    den = prod(factors)
    num = x*sum(s*prod(factors[:i]+factors[i+1:]) for i,s in enumerate(slopes))
    time_den = prod(1+j*j*x for j in range(1, 6))
    shell_den = spectrum*den
    shell_num = (spectrum+x*spectrum.derivative())*num
    contact = (x-left)*prod((x-v)**2 for v in nodes[1:3])
    cols = [divmod(x**i*shell_den, contact)[1] for i in range(5)]
    mat = arb_mat([[poly[i] for poly in cols] for i in range(5)])
    remainder = divmod(time_den*shell_num, contact)[1]
    fit = arb_poly(mat.solve(arb_mat([[remainder[i]] for i in range(5)])).entries())
    lead = (19*time_den(bound)*num(bound)/den(bound)-fit(bound))/contact(bound)
    numerator = fit+lead*contact
    quotient, residual = divmod(numerator*shell_den-time_den*shell_num, contact)
    assert all(v.contains(0) for v in residual)
    # 根为正的二次式确定 g 的变号处。重叠区间同时准备两个主控界。
    a0, a1, a2 = [p.derivative()(ref) for p in invariants]
    assert a0 < 0 < a1 and a2 > 0
    assert a0+73000*a1+73000**2*a2 < 0 < a0+74000*a1+74000**2*a2
    # x=(left+74000y)/(1+y) 将有界区间的判号化为正轴上的正系数。
    positive = sum(v*(left+74000*x)**i*(1+x)**(quotient.degree()-i)
                   for i, v in enumerate(quotient))
    assert all(v > 0 for v in positive)
    # 由低至高除以 bound-x，避免将 11 位数的端点反复乘入区间半径。
    lifted = numerator*den-19*time_den*num
    coefficients = []
    previous = arb(0)
    for i in range(8):
        previous = (lifted[i]+previous)/bound
        coefficients.append(previous)
    assert (lifted[8]+previous).contains(0)
    assert all(v > 0 for v in arb_poly(coefficients)(x+73000))
    sensitivity = arb(0)
    for j in range(1, 6):
        pole = -1/arb(j*j)
        others = prod(1+k*k*pole for k in range(1, 6) if k != j)
        sensitivity += abs(numerator(pole)/others)
    return margin, arb_poly(coefficients), time_den*den, contact, time_den, sensitivity, quotient, time_den*shell_den


def verify():
    ctx.prec = 1024
    bound = arb(117387924448)
    margin, _, _, _, _, sensitivity, _, _ = envelope(bound)
    assert margin+(1+sensitivity)*arb(2)**-81 < 0
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    scale = first*bound.sqrt()
    assert scale > arb("2899267716.58")
    print("尺度和导数的全域界", scale.str(20), margin.str(14))
    for cutoff, bits in ((10_000_000_000, 40), (50_000_000_000, 42), (100_000_000_000, 45)):
        error, _, _, _, _, sensitivity, _, _ = envelope(arb(cutoff))
        assert error+(1+sensitivity)*arb(2)**-bits < 0
        print(bits, "位绝对误差对应尺度下界", (first*arb(cutoff).sqrt()).str(15))
    print("若不使用首枚所需尺度", (second*bound.sqrt()).str(20))
    top = arb(2903616618)
    # 取第二枚零点在尺度上端的 x 上界。第一枚才可能越出该区。
    cap = ((top/second)**2).upper()
    deficit, high, den, contact, time_den, _, _, _ = envelope(cap)
    first_max = ((top/first)**2).upper()
    strip = (first_max-bound)/128
    rewards = []
    for j in range(128):
        interval = (bound+j*strip).union(bound+(j+1)*strip)
        rewards.append(((interval-cap)*high(interval)/den(interval)).upper())
    floor = -deficit/max(rewards)
    assert floor > arb("4.65271e-11")
    # 在首枚窗口中，核 B/Q 单调增加。其余节点对 B/Q 的贡献非负。
    times = model["shared"]["clock"]["clock"]["times"]()
    moments = neutral["contrast"]["moments_five"](times)
    expectation = sum(k*v for k,v in zip(contact.coeffs(),moments,strict=True))
    ceiling = expectation*time_den(bound)/contact(bound)
    assert ceiling < arb("4.66881e-11")
    print("小尺度解强制首枚的算术概率",floor.str(14),ceiling.str(14))


if __name__ == "__main__":
    verify()
