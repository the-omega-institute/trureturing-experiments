"""由六个参考区间的有理端点，求两个质量形状因子之比的精确极值。"""
from fractions import Fraction as Q
from math import prod

# ch39-two-zero-quark-window.py 中同一组 1 TeV 中心与半宽。
center = list(map(Q, ("0.00000615", "0.00311", "0.8616", "0.0000135", "0.000268", "0.01401")))
width = list(map(Q, ("0.00000014", "0.00005", "0.0043", "0.0000002", "0.000003", "0.00009")))
low = [value - error for value, error in zip(center, width, strict=True)]
high = [value + error for value, error in zip(center, width, strict=True)]


def sector_extremes(start):
    lower, upper = low[start:start+3], high[start:start+3]
    assert all(x > 0 for x in lower)
    # 单变量 s=y_i^2 满足 f^3 与 (s+A)^3/s 成正比，A 为其余两项平方和。
    # 0<s<t 时差商 s+t+3A-A^3/(st)，在二者均 ≤A/2 时为负，均 ≥A/2 时为正。
    for i in (0, 1):
        assert 2 * upper[i]**2 < sum(lower[j]**2 for j in range(3) if j != i)
    assert 2 * lower[2]**2 > upper[0]**2 + upper[1]**2
    return (upper[0], upper[1], lower[2]), (lower[0], lower[1], upper[2])


def shape_cube(masses):
    return sum(y**2 for y in masses)**3 / prod(masses)**2


up_min, up_max = sector_extremes(0)
down_min, down_max = sector_extremes(3)
lo = shape_cube(up_min) / shape_cube(down_max)
hi = shape_cube(up_max) / shape_cube(down_min)
assert Q("75.43570316720133")**3 < lo < Q("75.43570316720134")**3
assert Q("84.78930477373304")**3 < hi < Q("84.78930477373305")**3
print("两个端点均能达到，其立方恰为", lo, hi)
