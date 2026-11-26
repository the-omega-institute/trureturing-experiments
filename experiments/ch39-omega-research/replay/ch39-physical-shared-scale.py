"""三项物理窗口内，四个算术时间和参考味行列式共用尺度与正谱权。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


root = Path(__file__).parent
shared = run_path(str(root / "ch39-shared-flavor-clock-weights.py"))
physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
source = shared["source"]
indices = (8359, 1674503460, 1674725340, 6025261094, 6025461713)


def gradient(x, lambdas, slopes):
    return sum(s*x/(1 + z*x) for z, s in zip(lambdas, slopes, strict=True))


def upper_polynomial(limit, times, lambdas, slopes):
    # 在 x 的原坐标计算 dμ/Q 的矩。构造含右端点的三点求积。
    scaled = shared["clock"]["clock"]["cleared_moments"](times, 4)
    moments = [v/200**i for i, v in enumerate(scaled)]
    tilted = [limit*moments[i] - moments[i+1] for i in range(4)]
    gram = arb_mat([[tilted[i+j] for j in range(2)] for i in range(2)])
    constant, linear = gram.solve(arb_mat([[-tilted[2]], [-tilted[3]]])).entries()
    discriminant = linear**2 - 4*constant
    assert discriminant > 0
    left, middle = [(-linear + sign*discriminant.sqrt())/2 for sign in (-1, 1)]
    assert 0 < left < middle < limit
    nodes = [left, middle, limit]
    vandermonde = arb_mat([[x**i for x in nodes] for i in range(3)])
    weights = vandermonde.solve(arb_mat([[v] for v in moments[:3]])).entries()
    assert all(w > 0 for w in weights)
    for i in (3, 4):
        assert (sum(w*x**i for w, x in zip(weights, nodes)) - moments[i]).contains(0)

    coefficients = []
    for z, slope in zip(lambdas, slopes, strict=True):
        value = slope/((1 + z*left)**2 * (1 + z*middle)**2 * (1 + z*limit))
        for j in range(1, 5):
            value *= j*j - z
        coefficients.append(value)
    assert lambdas[0] < lambdas[1] < lambdas[2]
    assert coefficients[0] > 0 and coefficients[1] > 0 > coefficients[2]
    assert sum(coefficients) > 0
    score = arb(0)
    for x, w in zip(nodes, weights, strict=True):
        for j in range(1, 5):
            w *= 1 + j*j*x
        score += w*gradient(x, lambdas, slopes)
    return score, (left, middle), coefficients


def discrete_score(scale, times, heights, lambdas, slopes):
    weights = source["certify_weights"](times, scale, 4, heights)
    score = sum(w*gradient((scale/h)**2, lambdas, slopes)
                for w, h in zip(weights, heights, strict=True))
    return score, weights


def verify():
    ctx.prec, ctx.cap = 768, 3
    (e1, e2, e3), _, lambdas, slopes = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    ctx.prec = 224
    heights = [acb.zeta_zero(n).imag for n in indices]
    ctx.prec = 768
    lower = arb("41312660.31094738")
    limit = (lower/heights[0])**2
    upper, _, _ = upper_polynomial(limit, times, lambdas, slopes)
    assert upper < 0

    lo, hi = arb("41312675.40981774807"), arb("41312675.40981774808")
    assert discrete_score(lo, times, heights, lambdas, slopes)[0] < 0
    assert discrete_score(hi, times, heights, lambdas, slopes)[0] > 0
    scale = lo.union(hi)
    _, weights = discrete_score(scale, times, heights, lambdas, slopes)
    assert e1.degree() <= 1
    for h in heights:
        x = (scale/h)**2
        assert e2.derivative().derivative()(arb(0, 1))/x + e3.derivative().derivative()(arb(0, 1)) < 0

    flavor = physical["joint"]["flavor"]
    readout = flavor["reference_readout"]()
    probability = arb(0)
    for n, weight in zip(indices, weights, strict=True):
        _, alpha, curvature = physical["observables"](n)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curvature < arb("1.002")
        vub = readout(flavor["singlet"](n)[0])[2]
        probability += weight*vub**2
        print(n, "谱权、强耦合、曲率比", weight.str(12), alpha.str(12), curvature.str(12))
    reference = readout(flavor["singlet"](11661)[0])[2]
    deficit = 1 - probability.sqrt()/reference
    assert arb(".08000326") < deficit < arb(".08000327")
    print("共同尺度存在区间", scale.str(22))
    print("平均跃迁概率与等效幅度相对亏损", probability.str(14), deficit.str(12))


if __name__ == "__main__":
    verify()
