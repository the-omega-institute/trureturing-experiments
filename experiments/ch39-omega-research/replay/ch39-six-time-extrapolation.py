"""以六个算术时间给出后续每一项在连续径向模型中的精确上下界。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, arb_poly, ctx


def quadrature(endpoint, moments, denominator):
    shifted = [moments[j + 1] - endpoint * moments[j] for j in range(6)]
    matrix = arb_mat([[shifted[i + j] for j in range(3)] for i in range(3)])
    coeff = matrix.solve(arb_mat([[-a] for a in shifted[3:]])).entries()
    cubic = arb_poly([*coeff, 1])
    roots = cubic.complex_roots(tol=1e-55)
    assert len(roots) == 3 and all(z.imag.contains(0) for z in roots)
    inside = sorted((z.real for z in roots), key=float)
    assert 0 < inside[0] < inside[1] < inside[2] < 1
    nodes = [arb(endpoint), *inside]
    matrix = arb_mat([[v**j if j else 1 for v in nodes] for j in range(4)])
    weights = matrix.solve(arb_mat([[m] for m in moments[:4]])).entries()
    assert all(w > 0 for w in weights)
    for j in range(4, 7):
        assert sum(w * v**j for w, v in zip(weights, nodes)).overlaps(moments[j])
    return nodes, [w * denominator(v) for w, v in zip(weights, nodes)]


def readout(rule, k):
    nodes, probabilities = rule
    return sum(p * (1-y) / (1 + (k*k-1)*y) for y, p in zip(nodes, probabilities))


def verify(replay):
    radial = run_path(str(replay / "ch39-six-time-neutral-normalization.py"))
    ctx.prec = 384
    times = list(radial["neutral"]["model"]["shared"]["clock"]["clock"]["times"]()[:7])
    times[6] += arb("1.360900010e-8")
    denominator, moments = radial["compact_moments"](times)
    left, right = [quadrature(e, moments, denominator) for e in (0, 1)]
    offsets = [(1-y)/y for rule in (left, right) for y in rule[0][1:]]
    remote_mass = left[1][0]

    def width(k):
        s = k*k
        return remote_mass * prod(s-j*j for j in range(1, 7)) / prod(s+b for b in offsets)

    for k in (1, 2, 3, 4, 5, 6):
        for rule in (left, right):
            assert readout(rule, k).overlaps(times[k])
    for k in (7, 8, 10, 20, 50, 100, 1000, 10000, 100000):
        lo, hi = readout(right, k), readout(left, k)
        assert lo < hi and (hi-lo).overlaps(width(k))
        print(k, lo.str(18), hi.str(18), width(k).str(14))
    for limit in (".000001", ".001", ".01", ".1", ".25", ".35"):
        threshold = arb(limit)
        a, b = 6, 7
        while width(b) < threshold:
            a, b = b, 2*b
        while b-a > 1:
            mid = (a+b)//2
            if width(mid) < threshold:
                a = mid
            else:
                b = mid
        assert width(a) < threshold < width(b)
        print("首次超过宽度", limit, b)
    print("左端点质量", remote_mass.str(23))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
