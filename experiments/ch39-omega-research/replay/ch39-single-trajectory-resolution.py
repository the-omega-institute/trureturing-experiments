"""Resolve neighboring zeta heights on one fixed gauge trajectory with a floating scale."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

model = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))
q = model["q"]


def bounds():
    time = q("26.859482438")
    common = model["target"](10875, time)
    radius = q(".0002")
    local = [None] * 3
    full = q(".999").union(q("1.001"))

    def observe(i, inverse, rate):
        if i < 23*64:
            return
        log_rates = [r / x for r, x in zip(rate, inverse)]
        assert -2 < log_rates[0] < log_rates[1] < log_rates[2]
        if q(Fraction(i, 64)) > time + radius or q(Fraction(i+1, 64)) < time - radius:
            return
        for j, (deriv, ratio) in enumerate(zip(rate, model["ratios"])):
            value = deriv / ratio
            local[j] = value if local[j] is None else local[j].union(value)

    model["running"]["integrate"]([x * full for x in model["centers"]], 33*64, observe)
    rates = [x / (common + radius.union(-radius)*x) for x in local]
    # At the optimum the outer-channel signed log residuals are equal and opposite.
    slope = 2 + (rates[0] + rates[2])/2
    spread = (rates[2] - rates[0])/2
    for index, low, high in [(10874, "2.308", "2.311"), (10876, "1.134", "1.137")]:
        offset = (common / model["target"](index, time)).log()
        assert abs(offset) < radius * slope.lower()
        shift = -offset / slope
        error = abs(shift) * spread * 10**6
        assert q(low) < error < q(high)
        print(index, "optimal scale displacement", shift.str(10), "log-ppm", error.str(10))


if __name__ == "__main__":
    ctx.prec = 160
    bounds()
