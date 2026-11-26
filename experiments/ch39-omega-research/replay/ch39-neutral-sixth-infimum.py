"""保留真实辐边权重的一次中性壳读出对第六个算术时间的锐下界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx

root = Path(__file__).parent
sharp = run_path(str(root / "ch39-sharp-sixth-clock.py"))
neutral = run_path(str(root / "ch39-neutral-clock-response.py"))
model = neutral["model"]


def verify():
    ctx.prec, ctx.cap = 1024, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    moments = sharp["shared"]["moments_five"](times)
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    limit = sum(s/z for s, z in zip(slopes, poles, strict=True))
    spectrum = neutral["neutral_polynomial"]()
    spectrum = spectrum*(1/spectrum[0])

    def derivative(left):
        nodes, weights, tail = sharp["rule"](left, moments)
        return 19*tail*limit + sum(
            w*neutral["efficiency"](v, spectrum)*model["gradient"](v, poles, slopes)
            for w, v in zip(weights, nodes, strict=True))

    low, high = arb(0), arb("1e-5")
    for _ in range(140):
        mid = (low+high)/2
        score = derivative(mid)
        if score > 0:
            low = mid
        elif score < 0:
            high = mid
        else:
            break
    assert derivative(low) > 0 > derivative(high)
    nodes, weights, tail = sharp["rule"](low.union(high), moments)
    assert 0 < nodes[0] < nodes[1] < nodes[2] < 1 and all(w > 0 for w in [tail, *weights])
    assert (tail+sum(weights)).contains(1)
    readings = [sum(w*model["time_kernel"](v, j) for w, v in zip(weights, nodes, strict=True))
                for j in range(1, 7)]
    assert all(v.overlaps(t) for v, t in zip(readings[:5], times[1:6], strict=True))
    gap = readings[-1]-times[6]
    assert arb("1.33015563912468e-8") < gap < arb("1.33015563912470e-8")

    # H 的公共时间分母 Q 和味导数的分子、分母。
    x = arb_poly([0, 1])
    time_den = prod(1+j*j*x for j in range(1, 7))
    factors = [1+z*x for z in poles]
    flavor_den = prod(factors)
    flavor_num = x*sum(s*prod(factors[:i]+factors[i+1:]) for i, s in enumerate(slopes))
    den = spectrum*flavor_den
    num = (spectrum+x*spectrum.derivative())*flavor_num
    contact = prod((x-v)**2 for v in nodes)

    # 在重根多项式的商环内，六个余数系数确定 H 的六个自由系数。
    columns = [divmod(x**j*den, contact)[1] for j in range(6)]
    matrix = arb_mat([[column[i] for column in columns] for i in range(6)])
    remainder = divmod(time_den*(num-19*limit*den), contact)[1]
    fit = matrix.solve(arb_mat([[remainder[i]] for i in range(6)])).entries()
    numerator = 19*limit*time_den + arb_poly(fit)
    surplus, rem = divmod(numerator*den-time_den*num, contact)
    assert all(v.contains(0) for v in rem)
    # H 与单位辐边响应有相同的端值，原分子最高次精确消去。
    assert surplus[21].contains(0)
    assert all(surplus[j] > 0 for j in range(21))

    a, b, c = [p.derivative()(ref) for p in invariants]
    assert a < 0 < b and c > 0 and a+73000*b+73000**2*c < 0
    high_num = (numerator*flavor_den-19*time_den*flavor_num)(x+73000)
    assert high_num[9].contains(0) and all(high_num[j] > 0 for j in range(9))
    time_pole = -1/arb(36)
    c6 = numerator(time_pole)/prod(1+j*j*time_pole for j in range(1, 6))
    assert c6 > 0
    print("单位辐边接触公式的最小正节点", nodes[0].str(16))
    print("一次中性壳响应下第六项正偏差的下确界", gap.str(17))
    print("第六项本身的下确界", readings[-1].str(19))


if __name__ == "__main__":
    verify()
