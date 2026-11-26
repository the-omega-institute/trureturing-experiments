"""在实际联合零点的三塔角盒内，寻找半角坐标的最小公分母有理双钟。"""
from fractions import Fraction
from math import gcd
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def half_angle_box(angle, radius, scale):
    lo = ((angle - radius) / 2).tan()
    hi = ((angle + radius) / 2).tan()
    # 先向外包围，再对入围的有理数用原边界逐一严格复核。
    lower = int((lo * scale).lower().floor().unique_fmpz())
    upper = int((hi * scale).upper().ceil().unique_fmpz())
    return lower, upper, lo, hi


def numerators(interval, denominator, scale):
    lower, upper, _, _ = interval
    return range(-((-lower * denominator) // scale), upper * denominator // scale + 1)


def first_pair(boxes, scale, limit):
    for denominator in range(1, limit + 1):
        xs = numerators(boxes[0], denominator, scale)
        if not xs:
            continue
        for y in numerators(boxes[1], denominator, scale):
            for x in xs:
                for v, (_, _, low, high) in zip((x, y), boxes):
                    assert low < arb(v) / denominator < high
                return denominator, x, y
    raise AssertionError("给定范围内没有落入角盒的有理点")


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 256
    height = acb.zeta_zero(2_306_888_100_085).imag
    clocks = shell["geo"]["phases"](height)
    center = [clocks[k].arg() for k in (724, 804)]
    radius = arb("1e-9")
    scale = 2**96
    boxes = [half_angle_box(a, radius, scale) for a in center]
    d, x, y = first_pair(boxes, scale, 4_000_000)
    assert (d, x, y) == (1_044_988, -5_469, 72_843) and gcd(gcd(x, y), d) == 1
    print("最小半角公分母、两分子", d, x, y)
    for a, n in zip(center, (x, y)):
        re = Fraction(d*d - n*n, d*d + n*n)
        im = Fraction(2*d*n, d*d + n*n)
        assert re**2 + im**2 == 1
        arg = (acb(0, n) + d) / (d - acb(0, n))
        error = abs(arg.arg() - a)
        assert error < radius
        print("有理单位相位的实部、虚部、角误差", re, im, error.str(14))
    for bits in range(32):
        denominator = 2**bits
        choices = [numerators(interval, denominator, scale) for interval in boxes]
        if all(choices):
            nums = [values[0] for values in choices]
            for n, (_, _, lo, hi) in zip(nums, boxes):
                assert lo < arb(n) / denominator < hi
            assert (bits, nums) == (28, [-1_404_871, 18_711_836])
            print("二进制共同定点所需的最少小数位、两分子", bits, *nums)
            break
    else:
        raise AssertionError("未能选出二进制有理相位")


if __name__ == "__main__":
    verify()
