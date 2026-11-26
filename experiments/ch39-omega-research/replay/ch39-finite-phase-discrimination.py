"""以有理圆周参考和有限次受控回路，识别整个正零点谱上的指定编号。"""
from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    zeros = acb.zeta_zeros(11430,3)
    left,center,right = [(z-1)/z for z in zeros]
    assert arb.pi()/2 > left.arg() > center.arg() > right.arg() > 0
    normal = (left+right)/abs(left+right)
    gap = ((center-left)*normal.conjugate()).real
    low_share = ((center-right)*(left-right).conjugate()).real / abs(left-right)**2
    assert 0 < low_share < 1
    assert abs(center-low_share*left-(1-low_share)*right-gap*normal) < arb("1e-48")
    assert arb("2.22888797210e-17") < gap < arb("2.22888797212e-17")
    print("Single-return complex separation",gap.str(15),"weight at 11430",low_share.str(15))
    def rational_phase(code,bits):
        slope = arb(code)/(1 << bits)
        return acb(1,slope)/acb(1,-slope)
    one_step = rational_phase(49574543,40)
    assert right.arg()<one_step.arg()<left.arg()
    margin = min(((center-z)*one_step.conjugate()).real.lower() for z in (left,right))
    assert margin/gap > arb(".99998")
    print("40-bit single reference retains fraction", (margin/gap).str(12))
    phase = rational_phase(12690787451,48)
    assert right.arg() < phase.arg() < left.arg()
    offset = abs((center/phase).arg())
    minimum = min(abs((z/phase).arg()).lower() for z in (left,right))
    size = 1 << 32
    center_error = (size**2-1)*offset**2/12
    leakage = 1/(size*(minimum/2).sin())**2
    assert center_error < arb("1.55e-11") and leakage < arb(".009035")
    neighbor_distance = abs((center/right).arg())
    required = 2*(2*(arb(".99")*arb(".01")).sqrt()).acos()/neighbor_distance
    assert 559481889 < required < 559481890
    print("Controlled forward loops",size-1,"target loss below",center_error.str(12),
          "other acceptance below",leakage.str(12),"necessary queries",559481890)
    # 覆盖整个圆周阻带时的切比雪夫极值，以及达到相同整数次数的有理窗。
    nearest = abs((right/phase).arg())
    assert nearest < abs((left/phase).arg())
    optimal_stretch = 1/(nearest/2).cos()
    count = 1221962754
    critical = arb(10).acosh()/optimal_stretch.acosh()
    assert count-1 < critical < count
    stretch = 1 + arb(237690461445)/2**96
    assert 1 < stretch < optimal_stretch
    side = 1/(count*stretch.acosh()).cosh()**2
    loss = count**2*offset**2/4
    assert side < arb(".01") and loss < arb("3.76e-12")
    print("Continuous stopband minimum",count,"rational-window leakage",side.str(13),
          "target loss",loss.str(12))


if __name__ == "__main__":
    verify()
