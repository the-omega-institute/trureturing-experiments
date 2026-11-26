"""接近五项最小尺度的三窗零点概率必须先降后升。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

floor = run_path(str(Path(__file__).with_name("ch39-neutral-five-scale-floor.py")))


def verify():
    ctx.prec, ctx.cap = 768, 3
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    bound = arb(117387924448)
    scale_floor = first*bound.sqrt()
    scale_top = 2903600000
    # 沿用已证首点权重的核，算出第二枚最多能占多少剩余预算。
    cap = ((arb(2903616618)/second)**2).upper()
    _, _, _, contact, denominator, *_ = floor["envelope"](cap)
    times = floor["model"]["shared"]["clock"]["clock"]["times"]()
    moments = floor["neutral"]["contrast"]["moments_five"](times)
    expectation = sum(v*w for v, w in zip(contact.coeffs(), moments, strict=True))
    first_min = arb("4.65271e-11")
    remaining = expectation-first_min*contact(bound)/denominator(bound)
    second_x = bound*(first/second)**2
    second_max = remaining*denominator(second_x)/contact(second_x)
    assert 0 < second_max < arb("1.60921e-13")
    relative_efficiency = second_max/first_min*(second/first)**2
    assert relative_efficiency < arb(".006779")

    # 用四时间核的既有有限带证书给中间高度段的总概率下界。
    clock = floor["model"]["shared"]["clock"]["clock"]
    m = [v/200**j for j,v in enumerate(clock["cleared_moments"](times, 4))]
    left, right = arb(".00015"), arb(".015")
    residual = m[1]-left*m[0]-m[4]/right**3
    assert residual > 0
    mass = residual**2/m[2]
    assert mass > arb(".29045")
    lo, hi = arb(23672421776), arb(237077947240)
    assert lo < scale_floor/right.sqrt() and scale_top/left.sqrt() < hi
    def count_main(t):
        return t/(2*arb.pi())*((t/(2*arb.pi())).log()-1)+arb(7)/8
    count_upper = count_main(hi)-count_main(lo)+arb(".28")*(lo.log()+hi.log())
    third_min = mass/count_upper
    assert third_min > arb("3.6221e-13") and third_min/second_max > arb("2.25")

    fraction = (second/first).log()/(lo/first).log()
    assert 0 < fraction < 1 and third_min < first_min
    chord = first_min*(third_min/first_min)**fraction
    distortion = chord/second_max
    assert distortion > 259
    print("第二枚概率与逆平方效率比上界", second_max.str(13), relative_efficiency.str(12))
    print("较高处某枚概率下界及相对第二枚的倍率", third_min.str(13), (third_min/second_max).str(12))
    print("在对数高度上相对任意对数凹权重的必要校正倍率", distortion.str(13))


if __name__ == "__main__":
    verify()
