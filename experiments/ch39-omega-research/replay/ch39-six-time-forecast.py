"""六个时间读数给晚时刻预测留下的自由度，含六枚三窗零点的区间证书。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx

scale = 2_903_600_000
indices = (295873368963, 295873886415, 142269606336,
           142269938072, 349712678, 349952257)


def interpolatory_weights(nodes, denominator, moments):
    y = arb_poly([0, 1])
    probabilities = []
    for i, t in enumerate(nodes):
        others = nodes[:i] + nodes[i+1:]
        numerator = prod(y-v for v in others)
        expectation = sum(moments[k]*a for k, a in enumerate(numerator))
        probabilities.append(denominator(t)*expectation / prod(t-v for v in others))
    return probabilities


def left_radau(moments, denominator):
    matrix = arb_mat([[moments[i+j+1] for j in range(3)] for i in range(3)])
    p = arb_poly([*matrix.solve(arb_mat([[-v] for v in moments[4:7]])).entries(), 1])
    roots = p.complex_roots(tol=1e-80)
    assert len(roots) == 3 and all(z.imag.contains(0) for z in roots)
    nodes = [arb(0), *sorted([z.real for z in roots], key=float)]
    assert all(a < b for a, b in zip(nodes, nodes[1:])) and nodes[-1] < 1
    weights = interpolatory_weights(nodes, denominator, moments)
    assert all(w > 0 for w in weights)
    for j in range(4, 7):
        assert sum(w*t**j/denominator(t) for w, t in zip(weights,nodes)).overlaps(moments[j])
    return nodes, weights


def verify(root=Path(__file__).resolve().parent):
    clock = run_path(str(root / "ch39-arithmetic-height-bandwidth.py"))
    actual = run_path(str(root / "ch39-seven-actual-neutral-times.py"))
    compact = run_path(str(root / "ch39-six-time-neutral-normalization.py"))
    ctx.prec, ctx.cap = 512, 3
    times = clock["times"]()[:7]
    times[6] += arb("1.360900010e-8")
    denominator, moments = compact["compact_moments"](times)
    nodes, weights = left_radau(moments, denominator)
    assert arb(".39155545776518715") < weights[0] < arb(".39155545776518716")
    print("连续松弛左端点的权", weights[0].str(25), flush=True)
    constants = []
    for height, x, _ in actual["shell_nodes"](scale, indices, height_bits=176):
        constants.append(x/(1+x))
    ctx.prec = 512
    for t, lo, hi in zip(nodes[1:], (constants[1], constants[3], constants[5]),
                          (constants[0], constants[2], constants[4]), strict=True):
        assert lo < t < hi
    limiting = interpolatory_weights([arb(0), *constants], denominator, moments)
    assert all(w>0 for w in limiting)
    print("六枚实零点与极限列的七个权", *[v.str(24) for v in limiting], flush=True)
    assert arb(".39155545776046") < limiting[0] < arb(".39155545776049")

    # 这个区间检验同时覆盖所有足够高的有限候选。插值权重对这里的每个端点参数精确解矩方程。
    near = arb(0).union(arb("1e-28"))
    finite = interpolatory_weights([near, *constants], denominator, moments)
    assert min(w.lower() for w in finite) > 0
    assert finite[0] > arb(".39155545776046")
    print("远点 y 不超过 1e-28 时的七个权", *[v.str(16) for v in finite], flush=True)

    gap = weights[0] - limiting[0]
    assert 0 < gap < arb("4.715e-12")
    print("全体晚时刻两条界的宽度上限", gap.str(21), flush=True)
    for k in (7, 100, 10000, 1000000):
        exact_nodes = [arb(0), *constants]
        lower = sum(w*(1-y)/(1+(k*k-1)*y) for y,w in zip(exact_nodes, limiting))
        upper = sum(w*(1-y)/(1+(k*k-1)*y) for y,w in zip(nodes,weights))
        assert 0 < upper-lower < gap
        t = 1/arb(k*k-1)
        coefficient = (1+t)*prod(1-(j*j-1)*t for j in range(2, 7))
        remainder = sum(
            w/denominator(y)*coefficient*y/(t+y)
            *prod(((y-r)/(t+r))**2 for r in nodes[1:])
            for y,w in zip(constants, limiting[1:], strict=True))
        assert remainder.overlaps(upper-lower)
        print("时刻、实谱上确界的下上界", k, lower.str(21), upper.str(21), flush=True)


if __name__ == "__main__":
    verify()
