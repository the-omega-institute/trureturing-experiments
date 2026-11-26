"""对相邻真实零点的共同非负轮廓，验证支撑右端的连续下界。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    _, _, first, second = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["adjacent_fields"](root)
    slope = arb("5.62e-11")
    start = arb("1e-14")
    end = arb("8.39530e-9")
    pairs = []
    for (g, a, z), (other_g, other_a, h) in zip(first[1:], second[1:], strict=True):
        assert g == other_g and (a - other_a).contains(0)
        pairs.append((g, a, z, h))

    # 除以最轻的指数。零端附近其余 195 项均可用同一质量间隙控制。
    lead = min(range(len(pairs)), key=lambda j: float(pairs[j][3]))
    g, a, _, ref = pairs[lead]
    terms = [(s, gj, aj, m) for gj, aj, z, h in pairs for s, m in ((1, h), (-1, z))]
    del terms[2*lead]
    gap = min((m-ref).lower() for _, _, _, m in terms)
    assert gap > 0
    envelope = sum(abs(aj) + start*abs(gj)/(2*slope) for _, gj, aj, _ in terms)
    assert a - start*abs(g)/(2*slope) > envelope * (-gap/start).exp()

    def difference(x, z, h):
        # 抽出较轻端的指数，用 expm1 保存相邻两枚的质量差。
        if z > h:
            return (-(h-ref)/x).exp() * -((h-z)/x).expm1()
        if h > z:
            return (-(z-ref)/x).exp() * ((z-h)/x).expm1()
        return (z-h)/x * (-(h.union(z)-ref)/x).exp()

    def separated(x):
        return sum((aj+x*gj/(2*slope))*difference(x, z, h) for gj, aj, z, h in pairs) > 0

    count = 0
    left = start
    while left < end:
        right = min(2*left, end)
        pending = [(left, right)]
        while pending:
            lo, hi = pending.pop()
            if separated(lo.union(hi)):
                count += 1
            else:
                assert hi-lo > arb("1e-23"), (lo, hi)
                mid = (lo+hi)/2
                pending.extend(((lo, mid), (mid, hi)))
        left = right
    print("最轻质量间隙、正分离线性组合覆盖的闭区间数、右端",
          gap.str(12), count, end.str(10))


if __name__ == "__main__":
    verify()
