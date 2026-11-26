"""用同一归一正权认证算术时间前缀和参考味行列式的驻点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx


source = run_path(str(Path(__file__).with_name("ch39-flavor-scale-time-prefix.py")))
clock = source["clock"]
determinant = source["determinant"]


def tail_certificate(times):
    # 用五项读数限制 x >= 1000 的总概率，z = 200 x。
    moments = clock["clock"]["cleared_moments"](times, 5)
    gram = arb_mat([[moments[i + j + 1] for j in range(2)] for i in range(2)])
    coeff = gram.solve(arb_mat([[-moments[3]], [-moments[4]]]))
    polynomial = arb_poly([coeff[0, 0], coeff[1, 0], 1])
    threshold = arb(200000)
    discriminant = polynomial[1]**2 - 4 * polynomial[0]
    assert discriminant > 0
    roots = [(-polynomial[1] + sign * discriminant.sqrt()) / 2 for sign in (-1, 1)]
    assert 0 < roots[0] < roots[1] < threshold

    numerator = arb_poly([0, 1]) * polynomial**2
    denominator_at_threshold = arb(1)
    for j in range(1, 6):
        denominator_at_threshold *= 1 + j*j*threshold/200
    factor = denominator_at_threshold / numerator(threshold)
    mean = factor * sum(numerator[i]*moments[i] for i in range(6))
    sensitivity = factor * sum(abs(r) for r in clock["residues_of"](
        numerator, clock["partial_fractions"](5)))
    return mean, sensitivity


def flavor_data():
    invariant = determinant["invariant"]
    e1, e2, e3, _ = invariant["invariants"]()
    reference = invariant["flavor"]["singlet"](11661)[1]
    d1, d2, d3 = [p.derivative()(reference) for p in (e1, e2, e3)]
    mass_sum, mass_pairs, mass_product = [p(reference) for p in (e1, e2, e3)]
    lambdas = [arb(y)**2 for y in (".0000135", ".000268", ".01401")]
    slopes = [(d1*z*z - d2*z + d3)/(3*z*z - 2*mass_sum*z + mass_pairs) for z in lambdas]
    assert slopes[0] > 0 and slopes[1] > 0 and slopes[2] < 0
    return (e1, e2, e3), reference, lambdas, slopes


def verify():
    ctx.prec, ctx.cap = 768, 3
    (e1, e2, e3), reference, lambdas, slopes = flavor_data()
    d1, d2, d3 = [p.derivative()(reference) for p in (e1, e2, e3)]
    mass_sum, mass_pairs, mass_product = [p(reference) for p in (e1, e2, e3)]
    times = clock["clock"]["times"]()
    b = -d1
    linear, quadratic = d2 + d3 + b*mass_pairs, d3 + b*mass_product
    assert b > 0 and linear > 0 and quadratic > 0
    assert d1 + d2 + b*mass_sum + 1000*linear + 1000**2*quadratic < 0
    bound = b + sum(slope/z for slope, z in zip(slopes, lambdas) if slope > 0)
    tail, amplification = tail_certificate(times)
    error = (b*(1 - times[1]) - bound*tail)/(b + bound*amplification)
    assert error > arb("2.53705e-10")

    scale = arb.union(*determinant["scale_interval"]())
    indices = (1, 2529893, 2658737, 4035321, 9896608, 25215595)
    ctx.prec = 192
    heights = [acb.zeta_zero(n).imag for n in indices]
    ctx.prec = 768
    columns = []
    for height in heights:
        x = (scale/height)**2
        gradient = sum(slope*x/(1 + z*x) for slope, z in zip(slopes, lambdas))
        columns.append([arb(1), *[1/(1 + j*j*x) for j in range(1, 5)], gradient])
        assert e1.degree() <= 1 and e2.derivative().derivative()(arb(0, 1))/x + e3.derivative().derivative()(arb(0, 1)) < 0
    matrix = arb_mat(list(zip(*columns)))
    rhs = arb_mat([[v] for v in [*times[:5], 0]])
    weights = matrix.solve(rhs)
    assert all(w > 0 for w in weights.entries())
    assert all(v.contains(0) for v in (matrix*weights - rhs).entries())
    print("五项逐项最大误差下界", error.str(13))
    print("四项见证的六个正权", *(w.str(10) for w in weights.entries()))


if __name__ == "__main__":
    verify()
