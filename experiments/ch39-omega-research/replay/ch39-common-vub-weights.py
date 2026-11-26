"""四项时间、参考味真空与逐枚夸克跃迁概率的共同正权。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


model = run_path(str(Path(__file__).with_name("ch39-physical-shared-scale.py")))
indices = (8359, 1674503460, 1674725340, 6025261094, 6025461713, 1674353576, 6025296660)
scale = 41312680


def verify():
    ctx.prec, ctx.cap = 768, 3
    shared = model["shared"]
    (e1, e2, e3), reference_phase, lambdas, slopes = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    flavor = model["physical"]["joint"]["flavor"]
    readout = flavor["reference_readout"]()
    reference_vub = readout(flavor["singlet"](11661)[0])[2]
    assert (reference_vub-arb(".003732")).contains(0)
    # 七个已认证零点附近的二进有理中心。编号相邻两侧的真高度给出全谱间隙。
    ctx.prec = 384
    heights, centers, gaps, coordinates, probabilities = [], [], [], [], []
    for n in indices:
        before, height, after = [z.imag for z in acb.zeta_zeros(n-1, 3)]
        gap = min((height-before).lower(), (after-height).lower())
        assert gap > 0
        numerator = (height*2**192).floor().unique_fmpz()
        assert numerator is not None
        center = arb(numerator)/2**192
        assert abs(center-height) < arb(2)**-192
        heights.append(height)
        centers.append(center)
        gaps.append(gap)
        phase, coordinate = flavor["singlet"](n)
        coordinates.append(coordinate)
        probabilities.append(readout(phase)[2]**2)
        _, alpha, curvature = model["physical"]["observables"](n)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curvature < arb("1.002")
        ctx.prec = 384

    ctx.prec = 768
    xmax = (scale/heights[0])**2
    assert e1.degree() == 1 and e2.degree() == 2 and e2[2] < 0
    assert e2.derivative().derivative()(arb(0, 1))/xmax + e3.derivative().derivative()(arb(0, 1)) < 0
    columns = []
    for h, probability in zip(heights, probabilities, strict=True):
        x = (scale/h)**2
        columns.append([arb(1), *[1/(1+j*j*x) for j in range(1, 5)],
                        model["gradient"](x, lambdas, slopes), probability])
    matrix = arb_mat(list(zip(*columns, strict=True)))
    target = arb_mat([[y] for y in [*times[:5], 0, reference_vub**2]])
    weights = matrix.solve(target)
    assert all(w > 0 for w in weights.entries())
    assert all(v.contains(0) for v in (matrix*weights-target).entries())
    variance = sum(w*(c-reference_phase)**2 for w, c in zip(weights.entries(), coordinates, strict=True))
    assert variance > arb("1e-6")
    print("七项共用谱权", *(w.str(13) for w in weights.entries()))
    print("相位偏差平方均值及下型质量矩组合亏额", variance.str(12), (-e2[2]*variance).str(12))

    # 将七根谱线换成正有理函数。中央节点归一化以便界定响应矩阵的列扰动。
    grad_bound = sum(abs(s)/l for s, l in zip(slopes, lambdas, strict=True))
    feature_bound = max(arb(1), grad_bound.upper())
    inverse = matrix.inv()
    inverse_bound = max(sum(abs(inverse[i,j]).upper() for j in range(7)) for i in range(7))
    total = model["physical"]["mask"]["time_data"]()[0]
    width = arb(2)**-80
    column_error = arb(0)
    for h, center, gap in zip(heights, centers, gaps, strict=True):
        distance = gap-abs(center-h)
        assert distance > 0
        tail_bound = total*h**2*(width**2+(center-h)**2)/distance**2
        column_error += feature_bound*tail_bound
    eta = inverse_bound*column_error
    assert eta < 1
    deviation = eta/(1-eta)*max(w.upper() for w in weights.entries())
    assert min(w.lower() for w in weights.entries()) > deviation
    print("二进半宽与正规化响应系数最大改变量", width.str(10), deviation.str(12))


if __name__ == "__main__":
    verify()
