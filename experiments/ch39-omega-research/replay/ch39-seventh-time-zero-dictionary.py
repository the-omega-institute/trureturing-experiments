"""六项算术时间固定时，第七项在二十三枚字典上的锐界及全体三窗上的内外界。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_poly, ctx

scale = 2_903_600_000
lower_support = (8359, 143921293614265, 13219432687559, 295872524995,
                 295668865161, 142269938072, 142269606336)
upper_support = (143921298450585, 13219432687559, 295874398386, 295668865161,
                 142271930613, 142268163595, 8359)


def measure(poly, moments):
    return sum(coefficient*moments[k] for k, coefficient in enumerate(poly))


def certified_rule(support, coordinates, times, moments, denominator):
    y = arb_poly([0, 1])
    nodes = [coordinates[n] for n in support]
    probabilities = []
    for i, node in enumerate(nodes):
        others = nodes[:i]+nodes[i+1:]
        unnormalized = prod(y-z for z in others)
        probabilities.append(denominator(node)*measure(unnormalized, moments)/prod(node-z for z in others))
    assert all(w > 0 for w in probabilities)
    assert sum(probabilities).contains(1)
    for k in range(1, 7):
        reading = sum(w*(1-v)/(1+(k*k-1)*v) for w, v in zip(probabilities, nodes, strict=True))
        assert reading.overlaps(times[k])

    # (1+48y)(G-H)=q P。H 为六次插值多项式，G=(1-y)D/(1+48y)。
    factor = prod(y-v for v in nodes)
    numerator = (1-y)*denominator
    pole = -arb(1)/48
    q = numerator(pole)/prod(pole-v for v in nodes)
    assert q < 0
    interpolant, remainder = divmod(numerator-q*factor, 1+48*y)
    assert interpolant.degree() <= 6 and remainder[0].contains(0)
    result = measure(interpolant, moments)
    direct = sum(w*(1-v)/(1+48*v) for w, v in zip(probabilities, nodes, strict=True))
    assert direct.overlaps(result)
    return result, probabilities


def width_peak(coordinates, endpoint_rules):
    # 对 s=k²>36，插值余项的系数
    # -s*∏(s-j²)/∏(1+(s-1)v) 严格负，两个端点的支撑与权均不随 s 改变。
    delta = {}
    for sign, support, weights in endpoint_rules:
        for n, w in zip(support, weights, strict=True):
            delta[n] = delta.get(n, 0)+sign*w
    # z=s-36。以两支撑之并上的因子作公分母，该分母在 z≥0 严格正。
    shift = arb_poly([35, 1])
    factors = {n: 1+shift*coordinates[n] for n in delta}
    denominator = prod(factors.values())
    numerator = sum(w*(1-coordinates[n])*prod(f for i, f in factors.items() if i != n)
                    for n, w in delta.items())
    derivative = numerator.derivative()*denominator - numerator*denominator.derivative()
    assert derivative.degree() == 20
    assert all(derivative[i] > 0 for i in range(18))
    assert all(derivative[i] < 0 for i in range(18, 21))
    # 笛卡尔符号法则只允许一个正根，以下两端的异号保证它确实存在。
    a, b = arb("3489.58587069974"), arb("3489.58587069975")
    assert derivative(a*a-36) > 0 > derivative(b*b-36)
    adjacent = [numerator(k*k-36)/denominator(k*k-36) for k in (3489, 3490)]
    assert adjacent[0] < adjacent[1]
    assert arb("8.093669181490403e-6") < adjacent[1] < arb("8.093669181490406e-6")
    print("二十三枚上的最大区间宽度，其整数阶数为 3490", adjacent[1].str(17))


def endpoints(coordinates, times, moments, denominator):
    results = []
    rules = []
    for support, direction in ((lower_support, 1), (upper_support, -1)):
        nodes = [coordinates[n] for n in support]
        for n, v in coordinates.items():
            if n not in support:
                # q<0。所有其余列都满足严格对偶不等式，故该端点的权惟一。
                assert direction*prod(z-v for z in nodes) > 0
        result, probabilities = certified_rule(support, coordinates, times, moments, denominator)
        results.append(result)
        rules.append((-direction, support, probabilities))
    width_peak(coordinates, rules)
    lower, upper = results
    assert arb(".93341795164287272944314") < lower < arb(".93341795164287272944316")
    assert arb(".93341795164287312935473") < upper < arb(".93341795164287312935475")
    assert arb("3.9991158e-16") < upper-lower < arb("3.9991159e-16")
    print("二十三枚合格零点上的最小值、最大值及区间宽度",
          lower.str(23), upper.str(23), (upper-lower).str(12))
    return lower, upper


# 上端构造以 y=0 为极限节点，实际高度趋于无穷时七个权连续且仍严格为正。
spectrum_upper = (None, 295873886415, 295873368963, 142269938072, 142269606336,
                  350062491, 349501044)
spectrum_lower = (8359, 107140825305636, 143921298450585, 295870024485,
                  295869382715, 142269938072, 142269606336)


def spectrum_bounds(coordinates, times, moments, denominator, extrapolation):
    coordinates = {**coordinates, None: arb(0)}
    feasible_low, limiting_high = [
        certified_rule(support, coordinates, times, moments, denominator)[0]
        for support in (spectrum_lower, spectrum_upper)
    ]
    # 第 8359 枚是第一枚共同候选，故所有 y 位于这两个 Radau 端点之间。
    b = arb(".99999999999150664")
    assert coordinates[8359] < b < 1
    lower_rule, upper_rule = [extrapolation["quadrature"](e, moments, denominator) for e in (b, 0)]
    assert all(v < b for v in lower_rule[0][1:])
    low, high = [extrapolation["readout"](rule, 7) for rule in (lower_rule, upper_rule)]
    assert arb(".93341795164287258400990") < low < feasible_low < arb(".93341795164287258639827")
    assert arb(".93341795164501234318150") < limiting_high < high < arb(".93341795164501234976010")
    assert feasible_low-low < arb("2.389e-18") and high-limiting_high < arb("6.579e-18")
    print("全体合格零点上第七项下确界与上确界各自的内外界",
          *[v.str(23) for v in (low, feasible_low, limiting_high, high)])


def verify(root=Path(__file__).parent):
    pinned = run_path(str(root / "ch39-pinned-cp-zero-support.py"))
    radial = run_path(str(root / "ch39-six-time-neutral-normalization.py"))
    ctx.prec = 576
    times = list(radial["neutral"]["model"]["shared"]["clock"]["clock"]["times"]()[:7])
    times[6] += arb("1.360900010e-8")
    denominator, moments = radial["compact_moments"](times)
    indices = pinned["indices"]
    replay = pinned["joint"]["replay"]
    coordinates = {}
    for n, (_, x, _) in zip(indices, replay["shell_nodes"](scale, indices, height_bits=112), strict=True):
        coordinates[n] = x/(1+x)
    ctx.prec = 576
    endpoints(coordinates, times, moments, denominator)
    extra = sorted((set(spectrum_lower) | set(spectrum_upper))-set(indices)-{None})
    for n, (_, x, _) in zip(extra, replay["shell_nodes"](scale, extra, height_bits=112), strict=True):
        coordinates[n] = x/(1+x)
    ctx.prec = 576
    extrapolation = run_path(str(root / "ch39-six-time-extrapolation.py"))
    spectrum_bounds(coordinates, times, moments, denominator, extrapolation)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
