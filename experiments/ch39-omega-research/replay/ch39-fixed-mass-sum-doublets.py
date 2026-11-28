"""原五权、固定顶耦合与逐层线性质量和下，独立中性残差的最小分裂。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    ctx.prec = 192
    energy, portal = map(arb, ("1.176397e14", ".42995401872"))
    system = run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"](root, energy, portal)
    weights = system["weights_at"](portal)[0]
    data = system["original"]
    sites, p0, p1 = data["sites"], data["p0"], data["p1"]
    masses = [z + sign*portal*system["rho"] for z, sign in zip(system["masses"], system["signs"], strict=True)]
    assert len(masses) == 16 and min(masses) > 0
    lo, hi = arb(".3939"), arb(".3941")
    # K_x/K_∞ 是 (1+t/x)exp(-t/x) 对正权 dt/t² 的平均值。
    cap = max(m.upper() for m in masses)*(1+hi)**2 + p1/4
    assert weights[0] > 0 and weights[1] < 0 and weights[3] < 0
    assert weights[2] > 0 and weights[4] > 0 and sites[4] > sites[2]
    ratio = cap/sites[2]
    floor = (weights[2]+weights[4])*(1+ratio)*(-ratio).exp() + weights[1]+weights[3]
    assert floor > arb(".8")
    run_path(str(root/"ch39-neutral-residual-all-cutoffs.py"))["verify"]()
    light = arb("-.074256577194663365546").union(arb("-.074256577194663365544")) + arb(0, arb("4.907e-20"))
    gauge = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rule = gauge["angular_rule"](24, 57)
    log_ratio = (p1/p0).log()

    def layer_gains(split):
        gains = []
        for b in masses:
            gain = arb(0)
            for x, c in zip(sites, weights, strict=True):
                left, right = [gauge["angular_exposure"](1, b*(1+s)**2, x, p0, p1, log_ratio, rule)
                               for s in (-split, split)]
                gain += c*(left-right)/(48*arb.pi())
            gains.append(gain)
        return gains

    gains = layer_gains(hi)
    below, above = light + sum(layer_gains(lo)), light + sum(gains)
    smallest = min(g.lower() for g in gains)
    assert smallest > arb(".00035829056")
    assert above - smallest < arb("-.00032974240")
    print("最多十五层分裂的残差上界", (above-smallest).upper().str(12))
    assert below < arb("-2.14e-5") and above > arb("2.85e-5")
    assert (1 + lo)/(1 - lo) > arb("2.299")
    print("线性质量和守恒的全域斜率比较下界、两端中性残差", floor.lower().str(8), below.str(9), above.str(9))


if __name__ == "__main__":
    verify()
