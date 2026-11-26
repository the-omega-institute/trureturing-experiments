"""味读出效率倍率至多十九时，第六项算术时间的可检验偏差。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_poly, ctx

source = run_path(str(Path(__file__).with_name("ch39-forty-one-source-contrasts.py")))


def verify():
    ctx.prec, ctx.cap = 6144, 3
    times = source["shared"]["clock"]["clock"]["times"]()
    invariants, reference, poles, slopes = source["shared"]["flavor_data"]()
    _, _, _, roots = source["quadrature"](6, times)
    amplitudes = [s / z for s, z in zip(slopes, poles, strict=True)]
    limit = sum(amplitudes)
    residues = source["residues"](6, roots, poles, amplitudes)
    source["three_pole_sign"](residues)
    assert limit > 0 and 0 < poles[0] < poles[1] < poles[2] < 1
    # 驻点导数 g 在唯一正根以下为负。
    a, b, c = [p.derivative()(reference) for p in invariants]
    assert a < 0 < b and c > 0 and a + 73000*b + 73000**2*c < 0

    time_poles = [-1 / arb(j*j) for j in range(1, 7)]
    coefficients = []
    for i, pole in enumerate(time_poles):
        other = [y for j, y in enumerate(time_poles) if i != j]
        residue = prod(pole - r for r in roots) / prod(pole - y for y in other)
        bracket = 18*limit + sum(k/(1 + z*pole) for k, z in zip(residues, poles, strict=True))
        coefficients.append((i+1)**2 * residue * bracket)

    # 在 x >= 73000 上，清去正分母后直接核验 L - 19g 的符号。
    x = arb_poly([73000, 1])
    time_factors = [1 + j*j*x for j in range(1, 7)]
    flavor_factors = [1 + z*x for z in poles]
    time_den, flavor_den = prod(time_factors), prod(flavor_factors)

    def omit(factors, i):
        return prod(factors[:i]) * prod(factors[i+1:])

    excess = 19*limit*time_den*flavor_den
    excess += flavor_den*sum(v*omit(time_factors, i) for i, v in enumerate(coefficients))
    excess -= 19*time_den*x*sum(v*omit(flavor_factors, i) for i, v in enumerate(slopes))
    # 两个有理函数的右端值都是 19G，最高次项精确相消。
    assert excess[9].contains(0) and all(excess[j] > 0 for j in range(9))
    margin = -(19*limit + sum(v*t for v, t in zip(coefficients, times[1:7], strict=True)))
    assert margin > 0 and coefficients[-1] > 0
    sixth = margin/coefficients[-1]
    joint = margin/sum(abs(v) for v in coefficients)
    print("前五项精确时第六项高于参考值至少", sixth.str(17))
    print("六项最大绝对误差至少", joint.str(17))
    print("第六项参考值", times[6].str(18))
    assert sixth > arb("1.29166e-8") and joint > arb("1.22437e-11")
    # 若实际与参考都在同一二进格，按系数符号在对应格端点取上界。
    for bits in (35, 36):
        bins = 2**bits
        score = 19*limit
        for value, coefficient in zip(times[1:7], coefficients, strict=True):
            left = (bins*value).floor().unique_fmpz()
            assert left is not None and (coefficient > 0 or coefficient < 0)
            score += coefficient * (int(left) + int(coefficient > 0)) / bins
        if bits == 36:
            assert score < arb("-2.37054e-11")
        else:
            assert score > 0
        print(bits, "位同格时线性组合上界", score.str(12))



if __name__ == "__main__":
    verify()
