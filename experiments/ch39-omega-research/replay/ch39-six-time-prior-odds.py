"""六项给定算术时间限制每枚三窗零点的概率，以及相对逆平方先验的远近赔率。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, arb_poly, ctx

scale = 2_903_600_000


def radial_certificate(root):
    bands = run_path(str(root / "ch39-six-time-height-bands.py"))
    clock = run_path(str(root / "ch39-arithmetic-height-bandwidth.py"))
    ctx.prec, ctx.cap = 384, 3
    times = clock["times"]()[:7]
    times[6] += arb("1.360900010e-8")
    denominator, moments = bands["radial"]["compact_moments"](times)
    y = arb_poly([0, 1])

    def cutoff(height):
        return 1/(1+(arb(height)/scale)**2)

    c = cutoff(10**8)
    roots = [arb(v) for v in ("0.0000000145582674", ".0012125940914", ".0049020808845")]
    assert all(0 < r < c for r in roots)
    cubic = prod(y-r for r in roots)
    numerator = cubic**2 * (denominator(c)/cubic(c)**2)
    low_cap = sum(a*m for a,m in zip(numerator, moments, strict=True))
    assert arb("4.6768386650505e-11") < low_cap < arb("4.6768386650506e-11")

    def tail_edge(height, upper):
        cut = cutoff(height)
        nodes, weights = bands["cutoff_rule"](cut, moments, denominator)
        bands["verify_indicator"](nodes, denominator, upper)
        # 下界针对严格高于 height，上界包括截点本身。
        left = sum(w for t,w in zip(nodes, weights, strict=True) if t < cut)
        return left+weights[0] if upper else left

    band_floor = tail_edge(38*10**9, False)-tail_edge(45*10**9, True)
    tail_floor = tail_edge(13*10**12, False)
    assert band_floor > arb(".2343969") and tail_floor > arb(".2772896")
    print("低段概率上界、中段与远尾概率下界", low_cap.str(19), band_floor.str(19), tail_floor.str(19), flush=True)
    return low_cap, band_floor, tail_floor


def full_support_basis(root):
    pinned = run_path(str(root / "ch39-pinned-cp-zero-support.py"))
    rows, times, invariants, reference, probabilities, _ = pinned["joint"]["columns_at"](
        scale, pinned["lower_support"])
    ctx.prec = 512
    columns, target = pinned["problem"](rows, times, probabilities)
    basis = arb_mat(list(zip(*columns, strict=True)))
    weights = basis.solve(arb_mat([[v] for v in target])).entries()
    assert min(w.lower() for w in weights) > arb("4e-11")
    phase = arb(0,1)
    assert all(p(phase) > 0 for p in invariants)
    shift = arb_poly([reference, 1])
    displacement = phase-reference
    average = arb_poly([0])
    penalty = arb(0)
    for row,w in zip(rows, weights, strict=True):
        _, x, eta, _, _ = row
        det = 1 + x*invariants[0] + x**2*invariants[1] + x**3*invariants[2]
        ratio = (det*(1/det(reference)))(shift)
        average += w*eta*ratio
        if x>1000:
            divided = arb_poly(ratio.coeffs()[1:])(displacement)
            assert divided > 0
            maximum = max(arb(1), (1+displacement*divided).upper())
            penalty += w*eta*divided.lower()**2/(2*maximum**2)
    assert average[1].contains(0)
    gap = penalty-arb_poly(average.coeffs()[2:])(displacement)
    assert gap > arb("4e-13")
    print("十列基的最小概率、未归一对数行列式二次亏额的下界",
          min(w.lower() for w in weights).str(16), gap.lower().str(16), flush=True)


def verify(root):
    low_cap, band_floor, tail_floor = radial_certificate(root)
    short = run_path(str(root / "ch39-short-arc-born-weights.py"))
    data = short["records"]()
    assert len(set(short["indices"])) == len(data) == 15
    ctx.prec = 192
    heights = [row[0] for row in data]
    assert all(h < 10**8 for h in heights)
    assert all(a < b for a,b in zip(heights,heights[1:]))
    prior_low = sum(1/h**2 for h in heights)
    assert prior_low > arb("2.95816e-8")
    a, b, h = [arb(value) for value in (38*10**9, 45*10**9, 13*10**12)]
    two_pi = 2*arb.pi()

    def count_main(t):
        return t/two_pi*((t/two_pi).log()-1)+arb(7)/8

    count_cap = count_main(b)-count_main(a)+arb(".28")*(a.log()+b.log())
    assert count_cap < 25_189_334_981
    later_floor = band_floor/count_cap
    earlier_cap = low_cap/15
    assert later_floor > arb("2.98")*earlier_cap
    prior_tail = ((h/two_pi).log()+1)/(two_pi*h)+arb(".28")*(2*h.log()+arb(".5"))/h**2
    assert prior_tail < arb("3.59422e-13")
    odds = tail_floor*prior_low/(low_cap*prior_tail)
    assert odds > arb("4.87e14")
    print("中段全部零点数上界、后枚与前枚的概率比较系数",
          count_cap.str(17), (later_floor/earlier_cap).str(17), flush=True)
    print("未归一逆平方权的低段下界、远尾上界、远近赔率倍率下界",
          prior_low.str(17), prior_tail.str(17), odds.str(17), flush=True)
    full_support_basis(root)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
