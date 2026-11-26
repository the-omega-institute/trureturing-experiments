"""在原低能编号窗核验时间齐次的编号、链长混合更新。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
crossing = run_path(str(root / "ch39-gauge-crossing-short-chain.py"))
joint = crossing["joint"]
first, selected, last = 10815, 11431, 12092
length = 23226


def verify():
    ctx.prec = 96
    model = joint["flow"]["gauge"]["shortest"]["model"]
    offset = arb(model["fixed"]["OFFSET"].numerator) / model["fixed"]["OFFSET"].denominator
    time = arb("26.81676714125598")
    eps = arb("1.005e-7")
    radius = 2 * eps / arb(".1435")
    interval = time + arb(0, radius.upper())
    center = model["target"](selected, interval)
    error = eps + radius * (arb(".09") + 2 * center.upper())
    max_heat = arb("2e-17")
    common = arb(min(joint["leads"]["terminal_modes"](k)[0] for k in (1, 2, 3)))

    # 两个外邻点也参加检验，故编号寄存器本身可以延伸到窗外。
    heights = [z.imag for z in acb.zeta_zeros(first - 1, last - first + 3)]
    assert len(heights) == last - first + 3
    base = selected - first + 1
    factor = center * (heights[base]**2 + offset)
    reciprocals = [1 / (h*h + offset) for h in heights]
    angles = [2 * (1 / (2*h)).atan() for h in heights]
    margins = []
    for j in range(len(heights) - 1):
        assert heights[j] < heights[j + 1]
        step = factor * (reciprocals[j] - reciprocals[j + 1])
        distance_sum = factor * (abs(reciprocals[j] - reciprocals[base])
                                 + abs(reciprocals[j + 1] - reciprocals[base]))
        gauge_drop = 3 * common * step * (distance_sum - 2*error)
        phase_change = 6 * (angles[j] - angles[j + 1])
        prior_change = max_heat * ((heights[j + 1]**2 + arb(1)/4)
                                   / (heights[j]**2 + arb(1)/4)).log()
        margins.append(((gauge_drop - phase_change - prior_change) / max_heat).lower())
    assert min(margins) > 1000000
    print("相邻编号边数", len(margins), "最小向内能降的认证下界（热能单位）",
          min(margins).str(12), flush=True)

    # 在目标编号，原逐步界与长度配分和分别控制反向步和离开目标的单步。
    crossing["dynamics"]["verify"]()
    thermal = run_path(str(root / "ch39-thermal-plaquette-occupation.py"))
    thermal["verify"]()
    label_steps, remaining = 22000, 158000
    distance = max(selected - first, last - selected)
    assert 1 + label_steps < length
    label_exponent = 2 * (arb(label_steps)/4 - distance)**2 / label_steps
    growth_exponent = 2 * (arb(remaining)/4 - (length - 1))**2 / remaining
    total_steps = label_steps + remaining
    failure = ((-label_exponent).exp() + (-growth_exponent).exp()
               + total_steps * ((-arb(1000000)).exp() + (-arb(40000)).exp()
                                + (-arb(4000)).exp()))
    assert -failure.log() > 2100
    print("混合更新总步数", total_steps, "末时刻失败概率负对数下界", (-failure.log()).lower().str(12))


if __name__ == "__main__":
    verify()
