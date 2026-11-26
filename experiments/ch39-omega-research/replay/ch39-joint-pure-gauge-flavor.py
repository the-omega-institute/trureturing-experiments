"""在整个低能输入盒上，用两条守恒电流和四项纯态味读数识别零点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent


def gauge_band():
    ground = run_path(str(root / "ch39-floating-clock-ground.py"))
    mixing = ground["gauge"]["events"]["mixing"]
    model = mixing["model"]
    rates = ground["flow_bounds"]()
    low, high, logs = mixing["volume"]["input_bounds"]()
    anchor = arb("26.81676714125598")
    inverse, _, projected = model["evolve"](
        [a.union(b) for a, b in zip(low, high, strict=True)], logs, anchor, 4096
    )
    weak, color = model["ratios"][1:]
    assert color == 1
    mass_floor = 1 - 3 * arb("4e-5")
    tolerance = arb(2)**-21
    speed = (rates[2] - rates[1]).lower()
    assert speed > 0
    offset = (23 - anchor).union(33 - anchor)
    for _ in range(4):
        # 两支减法阈值皆非负。各支的谱权重归一，总质量有共同下界。
        z = (offset.upper().exp() / (42 * mass_floor))**2
        upper = min((z / 5).upper(), (1 + z / 4).log().upper()) / (4 * arb.pi()**2)
        radius = (abs(projected).upper() + 2*tolerance + upper*max(1/weak, 1/color)) / speed
        offset = arb(0, radius.upper())
    # 在最终尺度区间重新包围颜色阈值。
    z = (offset.upper().exp() / (42 * mass_floor))**2
    color_upper = z / (20 * arb.pi()**2)
    measured = inverse[2]/color + rates[2]*offset - arb(0).union(color_upper)/color
    measured += arb(0, tolerance)
    assert measured > 0
    fixed = model["fixed"]
    factor = (fixed["geo"]["BOOK_PLANCK"]/1000)**2 * 3*arb(5).sqrt()/(4*arb.pi())
    heights = (factor * (-2*(anchor + offset)).exp() / measured - model["q"](fixed["OFFSET"])).sqrt()
    counts = [edge.zeta_nzeros().unique_fmpz() for edge in (heights.lower(), heights.upper())]
    assert all(count is not None for count in counts)
    first, last = int(counts[0]) + 1, int(counts[1])
    assert (first, last) == (8952, 14578)
    assert abs(offset).upper() < arb(".208116")
    # 严格的端点比较也处理闭容差及可能贴边的零点。
    assert acb.zeta_zero(first-1).imag < heights.lower() < acb.zeta_zero(first).imag
    assert acb.zeta_zero(last).imag < heights.upper() < acb.zeta_zero(last+1).imag
    print("Gauge scale offset below", abs(offset).upper().str(10), "outer labels", first, last, flush=True)
    return first, last


def verify():
    ctx.prec = 160
    first, last = gauge_band()
    evaluate = run_path(str(root / "ch39-high-pure-flavor-counterexamples.py"))["flavor_test"]()
    accepted = []
    runner = None
    for start in range(first, last+1, 256):
        count = min(256, last-start+1)
        zeros = acb.zeta_zeros(start, count)
        assert len(zeros) == count
        for n, zero in enumerate(zeros, start):
            errors = evaluate(zero)
            if all(value < 1 for value in errors):
                accepted.append(n)
            else:
                lower = max(value.lower() for value in errors)
                assert lower > 1
                if runner is None or lower < runner[0]:
                    runner = (lower, n)
    assert accepted == [11431]
    assert runner[0] > arb("1.79846")
    print("All four pass at", accepted, "minimum competing lower bound", runner[0].str(12), "at", runner[1])


if __name__ == "__main__":
    verify()
