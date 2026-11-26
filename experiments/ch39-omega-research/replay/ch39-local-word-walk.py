"""有限对称字游走的末圈份额、首达步数及最大停步概率。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    geometry = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 224
    theta = geometry["unit"]((108469, 109)).arg()
    assert 0 < 106 * theta < arb.pi()
    cosines = [(q * theta).cos() for q in range(1, 106)]
    assert all(a > b for a, b in zip(cosines, cosines[1:]))
    assert cosines[-1] < 0
    target = arb(9) / 125

    def shares(stay):
        factors = [(stay + (1 - stay) * c)**2 for c in cosines]
        powers = [arb(1) for _ in factors]
        values = []
        for _ in range(107):
            values.append(powers[-1] / sum(powers))
            powers = [a*b for a, b in zip(powers, factors)]
        return values

    fair = shares(arb(0))
    assert all(x < target for x in fair[:17])
    assert all(x > target for x in fair[17:])
    print("无停步的第 16、17、105、106 步", *[fair[k].str(16) for k in (16, 17, 105, 106)])

    low = arb("0.004824685981972")
    high = arb("0.004824685981974")
    left, right = shares(low), shares(high)
    assert left[38] > target > right[38]
    assert all(val < target for k, val in enumerate(left) if k != 38)
    print("极限停步概率的两端、第 38 步份额", low.str(16), high.str(16),
          left[38].str(17), right[38].str(17))

    for denom in (213, 867):
        values = shares(arb(1) / denom)
        passing = [k for k, v in enumerate(values) if v > target]
        assert all(v < target or v > target for v in values)
        print("停步概率分母、过门槛的首末步数、第 106 步份额",
              denom, passing[0], passing[-1], values[-1].str(16))
        if denom == 867:
            assert values[-1] > arb(".14243")

    # 有限时域的最优停止。每层只保留从中心在剩余预算内可达的站。
    payoff = [(5 * theta * j).cos() for j in range(-106, 107)]
    values = payoff
    for step in range(1, 107):
        values = [stop.min((left + right) / 2)
                  for stop, left, right in zip(payoff[step:-step], values[:-2], values[2:], strict=True)]
    lower, = values
    assert arb(".25496818870") < lower < arb(".25496818871")
    assert lower**2 / 105 > arb(".000619")
    print("106 步以内任选停时的第五矩实部最小值", lower.str(17))


if __name__ == "__main__":
    verify()
