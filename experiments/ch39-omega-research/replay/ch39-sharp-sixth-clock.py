"""前五项精确且味效率倍率至多十九时，第六项的锐下界。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

shared = run_path(str(Path(__file__).with_name("ch39-five-time-source-contrast.py")))
model = shared["model"]


def rule(left, moments):
    # 尾节点只贡献第五个清分母矩。其余矩由含指定左端的三点公式给出。
    shifted = [b - left*a for a, b in zip(moments[:4], moments[1:5], strict=True)]
    gram = arb_mat([[shifted[0], shifted[1]], [shifted[1], shifted[2]]])
    p, q = gram.solve(arb_mat([[-shifted[2]], [-shifted[3]]])).entries()
    gap = (q*q - 4*p).sqrt()
    nodes = [left, (-q-gap)/2, (-q+gap)/2]
    matrix = arb_mat([[x**j for x in nodes] for j in range(3)])
    tilted = matrix.solve(arb_mat([[v] for v in moments[:3]])).entries()
    weights = [w*prod(1+j*j*x for j in range(1, 6)) for w, x in zip(tilted, nodes, strict=True)]
    tail = 14400*(moments[5] - sum(w*x**5 for w, x in zip(tilted, nodes, strict=True)))
    return nodes, weights, tail


def verify():
    ctx.prec, ctx.cap = 768, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    moments = shared["moments_five"](times)
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    amplitudes = [s/z for s, z in zip(slopes, poles, strict=True)]
    limit = sum(amplitudes)

    def derivative(left):
        nodes, weights, tail = rule(left, moments)
        return sum(w*model["gradient"](x, poles, slopes) for x, w in zip(nodes, weights, strict=True)) + 19*tail*limit

    low, high = arb(0), arb("1e-5")
    assert derivative(low) > 0 > derivative(high)
    for _ in range(135):
        mid = (low+high)/2
        value = derivative(mid)
        if value > 0:
            low = mid
        elif value < 0:
            high = mid
        else:
            break
    assert derivative(low) > 0 > derivative(high)
    nodes, weights, tail = rule(low.union(high), moments)
    assert 0 < nodes[0] < nodes[1] < nodes[2] < 1 and all(v > 0 for v in [tail, *weights])
    assert (tail + sum(weights)).contains(1)
    def reading(j):
        return sum(w*model["time_kernel"](x, j) for x, w in zip(nodes, weights, strict=True))
    assert all(reading(j).overlaps(times[j]) for j in range(1, 6))
    lower_gap = reading(6) - times[6]
    assert arb("1.32318575047378e-8") < lower_gap < arb("1.32318575047379e-8")

    # 接触根成双。低 x 段的非负余项由三极点次序给出。
    a, b, c = [p.derivative()(ref) for p in invariants]
    assert a < 0 < b and c > 0 and a + 73000*b + 73000**2*c < 0
    assert limit > 0 and 0 < poles[0] < poles[1] < poles[2] < 1
    roots = [x for x in nodes for _ in range(2)]
    residues = [v*prod(1-z/(j*j) for j in range(1, 7))/prod(1+z*x for x in roots)
                for v, z in zip(amplitudes, poles, strict=True)]
    assert amplitudes[0] > 0 and amplitudes[1] > 0 > amplitudes[2]
    assert residues[0] > 0 and residues[1] > 0 > residues[2] and sum(residues) > 0
    # 在高 x 段，G-g >= G/(1+λ₃x)，1-S <= (Σr + Σ1/j²)/x。
    assert 73000/(1+73000*poles[2]) > sum(roots) + sum(1/arb(j*j) for j in range(1, 7))
    pole6 = -1/arb(36)
    coefficient6 = 36*prod(pole6-r for r in roots)/prod(pole6+1/arb(j*j) for j in range(1, 6))
    bracket = 18*limit + sum(k/(1+z*pole6) for k, z in zip(residues, poles, strict=True))
    assert coefficient6*bracket > 0
    print("最小内部节点", nodes[0].str(20))
    print("前五项精确时第六项正偏差下确界", lower_gap.str(20))
    print("第六项的下确界", (times[6]+lower_gap).str(20))


if __name__ == "__main__":
    verify()
