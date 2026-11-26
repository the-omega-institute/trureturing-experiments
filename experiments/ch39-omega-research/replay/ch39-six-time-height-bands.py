"""六项给定算术时间在三个实际零点高度区间内各要求多少谱权。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx

radial = run_path(str(Path(__file__).with_name("ch39-six-time-neutral-normalization.py")))
y = arb_poly([0, 1])
scale = 2_903_600_000


def cutoff_rule(cutoff, moments, denominator):
    fixed = [cutoff, arb(0), arb(1)]
    factor = prod(y-v for v in fixed)
    modified = [sum(moments[i+j]*a for j, a in enumerate(factor)) for i in range(4)]
    pair = arb_mat([[modified[i+j] for j in range(2)] for i in range(2)])
    coeff = pair.solve(arb_mat([[-v] for v in modified[2:]])).entries()
    roots = arb_poly([*coeff, 1]).complex_roots(tol=1e-52)
    assert len(roots) == 2 and all(v.imag.contains(0) and 0 < v.real < 1 for v in roots)
    nodes = [*fixed, *sorted((v.real for v in roots), key=float)]
    assert nodes[3] < nodes[4] and all(v < cutoff or v > cutoff for v in nodes[3:])
    matrix = arb_mat([[v**k if k else 1 for v in nodes] for k in range(5)])
    weights = matrix.solve(arb_mat([[v] for v in moments[:5]])).entries()
    assert all(w > 0 for w in weights)
    for k in (5, 6):
        assert sum(w*v**k for w, v in zip(weights, nodes)).overlaps(moments[k])
    return nodes, [w*denominator(v) for w, v in zip(weights, nodes)]


def positive_quotient(poly, factors, left, right):
    # 已知插值零点恰可除尽；保留商的外包围。每个除去的因子在检验侧非负。
    quotient, remainder = divmod(poly, prod(factors))
    assert all(v.contains(0) for v in remainder)
    pending = [(left, right)]
    while pending:
        a, b = pending.pop()
        mid = (a+b)/2
        if quotient(arb_poly([mid, 1]))(arb(0, b-mid)) > 0:
            continue
        assert b-a > arb("1e-12")
        pending.extend(((a, mid), (mid, b)))


def verify_indicator(nodes, denominator, upper):
    cutoff = nodes[0]
    matrix, values = [], []
    for i, v in enumerate(nodes):
        on_left = v < cutoff or (i == 0 and upper)
        matrix.append([v**k if k else 1 for k in range(7)])
        values.append([denominator(v) if on_left else 0])
        if i >= 3:
            matrix.append([0, 1, *[k*v**(k-1) for k in range(2, 7)]])
            values.append([denominator.derivative()(v) if on_left else 0])
    polynomial = arb_poly(arb_mat(matrix).solve(arb_mat(values)).entries())
    left = [y, *[(y-v)**2 for v in nodes[3:] if v < cutoff]]
    right = [1-y, *[(y-v)**2 for v in nodes[3:] if v > cutoff]]
    if upper:
        left.append(cutoff-y)
        positive_quotient(polynomial-denominator, left, arb(0), cutoff)
        positive_quotient(polynomial, right, cutoff, arb(1))
    else:
        right.append(y-cutoff)
        positive_quotient(denominator-polynomial, left, arb(0), cutoff)
        positive_quotient(-polynomial, right, cutoff, arb(1))


def verify():
    ctx.prec = 384
    neutral = radial["neutral"]
    times = list(neutral["model"]["shared"]["clock"]["clock"]["times"]()[:7])
    times[6] += arb("1.360900010e-8")
    denominator, moments = radial["compact_moments"](times)
    heights = ("3.8e10", "4.5e10", "7e10", "1e11", "1e12", "2.4e13")
    tails = {}
    for h in heights:
        cutoff = 1/(1+(arb(h)/scale)**2)
        nodes, weights = cutoff_rule(cutoff, moments, denominator)
        verify_indicator(nodes, denominator, False)
        verify_indicator(nodes, denominator, True)
        below = sum(w for v, w in zip(nodes, weights) if v < cutoff)
        tails[h] = below.union(below+weights[0])
    bands = [tails[heights[i]]-tails[heights[i+1]] for i in (0, 2)] + [tails[heights[4]]]
    assert tails[heights[5]] > arb(".0021052")
    assert sum(bands) > arb(".99915")
    polynomial = neutral["neutral_polynomial"]()
    normal_lo = arb("1.0020816770805")
    normal_hi = arb("1.0020816771057")
    change = 18*arb(23)/25*arb(2)**-20/arb(8462)**2
    flavor = []
    for i, mass in enumerate(bands):
        xmax = (scale/arb(heights[2*i]))**2
        xmin = (scale/arb(heights[2*i+1]))**2 if i < 2 else arb(0)
        response = neutral["efficiency"](xmin, polynomial).union(neutral["efficiency"](xmax, polynomial)+change)
        flavor.append(mass*response/normal_lo.union(normal_hi))
    assert sum(flavor) > arb(".998")
    for title, ranges in (("算术概率", bands), ("味概率", flavor)):
        print(title, *[(v.lower().str(12), v.upper().str(12)) for v in ranges])
    print("三段合计的算术、味概率下界与最高远尾下界",
          sum(bands).lower().str(12), sum(flavor).lower().str(12), tails[heights[5]].lower().str(12))


if __name__ == "__main__":
    verify()
