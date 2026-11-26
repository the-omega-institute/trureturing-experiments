"""实际十枚正权表的 19 项最小时间辨识数及有限精度总变差证书。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx

scale = 2_903_600_000
between = (349501044, 142268184716, 142269938072, 142270888235, 295869382715,
           295873352222, 295873492883, 295874150547, 107140825305636, 143921294494500)


def denominator(rates, x):
    return prod(1+s*x for s in rates)


def partial_cost(roots, rates, multiplier):
    # 正规化精确已知，独立时间误差只作用于极点对应的系数。
    result = arb(0)
    for s in rates:
        pole = -arb(1)/s
        residue = multiplier*prod(pole-r for r in roots)/prod(1+t*pole for t in rates if t != s)
        result += abs(residue)
    return result


def stability(indices):
    ctx.prec = 112
    neighbors = {n: [z.imag for z in acb.zeta_zeros(n-1, 3)] for n in indices}
    ctx.prec = 640
    triples = [[(scale/h)**2 for h in neighbors[n]] for n in indices]
    assert all(a > b > c > 0 for a, b, c in triples)
    nodes = [row[1] for row in triples]
    assert all(a > b for a, b in zip(nodes, nodes[1:]))
    roots = [nodes[0], *[v for z in nodes[1:] for v in (z, z)]]
    rates = [j*j for j in range(1, 20)]

    def exposing(x):
        return (nodes[0]-x)*prod((x-z)**2 for z in nodes[1:])/denominator(rates, x)

    # 每对根之间，log(Q(exp(t))) 对 t 严格凹，故网格最小值只需考察邻点及 0。
    candidates = [arb(0), triples[0][2], *[v for a, _, c in triples[1:] for v in (a, c)]]
    gaps = [exposing(v) for v in candidates]
    gap = min(v.lower() for v in gaps)
    assert gap > 0
    exposure_cost = partial_cost(roots, rates, arb(-1))
    leakage = exposure_cost/gap
    # 十个有理 Lagrange 函数在支撑上分别为单点指示函数。
    # 每段上 log(|Lᵢ|/Q) 关于 log x 严格凸，所以它在网格上的最大值也只在这些端点中取。
    cardinal_cost = arb(0)
    ratios = arb(0)
    for i, z in enumerate(nodes):
        others = nodes[:i]+nodes[i+1:]
        multiplier = denominator(rates[:9], z)/prod(z-a for a in others)
        cardinal_cost += partial_cost(others, rates[:9], multiplier)
        values = [abs(multiplier*prod(v-a for a in others))/denominator(rates[:9], v)
                  for v in candidates]
        ratios += max((value/q).upper() for value, q in zip(values, gaps, strict=True))
    tv = (leakage+cardinal_cost+ratios*exposure_cost)/2
    exponent = int((tv.log()/arb(2).log()).ceil().unique_fmpz())
    assert exponent == 232 and tv < arb(2)**232
    print("名单外质量与总变差的误差放大上界（以 2 为底的对数）",
          (leakage.log()/arb(2).log()).str(14), (tv.log()/arb(2).log()).str(14), exponent)
    return exponent


def verify(root):
    pinned = run_path(str(root / "ch39-pinned-cp-zero-support.py"))
    joint = pinned["joint"]
    base = sorted(pinned["previous"]["indices"][:10])
    assert base[0] == 8359
    assert all(a < i < b for a, i, b in zip(base[:-1], between[:-1], base[1:], strict=True))
    assert base[-1] < between[-1]
    # 既有十枚的严格正权见 ch39-eleven-neutral-cp-range.py 中 pinned_range。
    # 十九项上界对所有这十枚上的正权都成立。只需另验交错处的十枚确实合格。
    list(joint["replay"]["shell_nodes"](scale, between, height_bits=112))
    stability(base)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
