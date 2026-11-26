"""用一只弧差钟在连续两枚实际零点间的相位增量，识别两段有理球心角。"""
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def ball(rational):
    return arb(rational.numerator)/rational.denominator


def verify():
    root = Path(__file__).resolve().parent
    graph = run_path(str(root/"ch39-arc-relative-determinant.py"))
    tomography = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    faces = tomography["enumerate_faces"]()
    ctx.prec = 192
    radius_squared = graph["R"]
    assert radius_squared == 867
    # 先由原距离推出各有理余弦。
    reference = 1-Fraction(244, 2*radius_squared)
    assert reference == Fraction(745, 867)
    first = 2_306_888_100_085
    gaps = []
    for index, max_gap in ((2, 4), (first, 2)):
        h0, h1 = [v.imag for v in acb.zeta_zeros(index, 2)]
        gap = h1-h0
        assert 0 < gap < max_gap
        gaps.append((index, gap, graph["phases"](h0), graph["phases"](h1)))
    for length, numerator in ((724, 505), (804, 465)):
        tip = 1-Fraction(length, 2*radius_squared)
        assert tip == Fraction(numerator, 867)
        product = reference * tip
        squared_radical = (1-reference**2)*(1-tip**2)
        assert any(isqrt(n)**2 != n for n in (squared_radical.numerator, squared_radical.denominator))
        cosine = ball(product) + ball(squared_radical).sqrt()
        constant = reference**2 + tip**2 - 1
        assert (cosine*cosine - 2*ball(product)*cosine + ball(constant)).contains(0)
        omega = ball(tip).acos() - ball(reference).acos()
        assert (omega.cos()-cosine).contains(0)
        assert reference.denominator // gcd(reference.denominator, 3) == 17**2
        face = tomography["axis_x" if length == 724 else "axis_y"][0]
        frequency, label = faces[face]
        assert frequency == ((1, 0) if length == 724 else (0, 1))
        bare = tomography["unit"](label)
        for index, gap, phases0, phases1 in gaps:
            observed = phases1[length]/phases0[length]
            assert (observed - acb(0, gap*omega).exp()).contains(0)
            assert 0 < gap*omega < arb.pi()
            if index == first:
                w0, w1 = bare*phases0[length], bare*phases1[length]
                cross_reading = (w1*w0.conjugate()).real
                recovered = cross_reading.acos()/gap
                assert (recovered-omega).contains(0)
                print("原壳单钟三角面、跨高度 Wilson 实部", face, cross_reading.str(19))
            print("钟弦、起始编号、高度差、解回的弧差频率",
                  length,index,gap.str(15),(observed.arg()/gap).str(16))


if __name__ == "__main__":
    verify()
