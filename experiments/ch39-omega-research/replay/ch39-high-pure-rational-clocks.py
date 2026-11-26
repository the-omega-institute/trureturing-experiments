"""用有限有理钟码、矩阵与有理数比较重放三个纯态的 CKM 检验。"""
from fractions import Fraction
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import fmpq as Q

finite = run_path(str(Path(__file__).with_name("ch39-rational-local-flavor.py")))
clock_codes = {
    11431: ((-1, -496481263766), (-1, 489993214324)),
    79545: ((-1, -470725993259), (-1, 517116220041)),
    97370: ((-1, -492842502141), (-1, 494817723896)),
}


def tangent_bounds(angle):
    # 在两个正端点小于 1.2 弧度时，以下相邻交错部分和夹住正弦与余弦。
    def partial(parity, count):
        return sum((-1)**k * angle**(2*k + parity) / factorial(2*k + parity)
                   for k in range(count))
    sin_low, sin_high = partial(1, 6), partial(1, 7)
    cos_low, cos_high = partial(0, 6), partial(0, 7)
    assert 0 < sin_low < sin_high and 0 < cos_low < cos_high
    return sin_low / cos_high, sin_high / cos_low


def verify():
    centers = [Q(22501, 100000), Q(4183, 100000), Q(3732, 1000000)]
    radius = Q(55, 1000000)
    _, low = tangent_bounds(Fraction(1141, 1000))
    high, _ = tangent_bounds(Fraction(1153, 1000))
    for n, codes in clock_codes.items():
        squares, jcp, cos_term = finite["measurements"](finite["flavor_operation"](codes))
        assert all((center-radius)**2 < value < (center+radius)**2
                   for center, value in zip(centers, squares, strict=True))
        assert jcp > 0 and cos_term > 0
        tangent = jcp / ((1 - squares[2]) * cos_term)
        assert Q(low.numerator, low.denominator) < tangent < Q(high.numerator, high.denominator)
        print("Rational clock codes pass all four flavor tolerances at", n)


if __name__ == "__main__":
    verify()
