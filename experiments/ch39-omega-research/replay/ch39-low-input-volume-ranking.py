"""Bound the 11488/10875 tube-volume ratio over the full low-energy Yukawa box."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
gauge = run_path(str(root / "ch39-floating-scale-winner.py"))
spectral = run_path(str(root / "ch39-hardy-height-selection.py"))
q = gauge["q"]
loop = gauge["loop"]
running = gauge["running"]
ratios = gauge["ratios"]


def gauge_derivative_defect(low, logs):
    # Bounds on row sums of D_{X(0)} X and D_{X(0)} Y with fixed initial Y.
    x_norm, y_norm = arb(1), arb(0)

    def observe(step, gauge_squares, traces):
        nonlocal x_norm, y_norm
        assert all(0 < g < 2 for g in gauge_squares)
        assert all(0 < t < 1 for t in traces)
        fourth = [g**2 for g in gauge_squares]

        def row_bound(matrix, factors):
            return max(sum(abs(k)*v for k, v in zip(row, factors)).upper()
                       for row in matrix)

        xx = 2 / loop**2 * row_bound(running["B"], fourth)
        xy = 2 / loop**2 * row_bound(running["D"], traces)
        yx = 2 / loop * row_bound(running["C"], fourth)
        total = 3*(traces[0] + traces[1]) + traces[2]
        yy = (3*sum(traces) + 2*total) / loop
        norm = max((xx + xy).upper(), (yx + yy).upper())
        # Increment of either variational block within this step (Gronwall).
        change = step*norm*(step*norm).exp()*max(x_norm, y_norm)
        x = x_norm + step*(xx*x_norm + xy*y_norm + (xx + xy)*change)
        y = y_norm + step*(yx*x_norm + yy*y_norm + (yx + yy)*change)
        x_norm, y_norm = x.upper(), y.upper()

        # Both endpoint exclusions, and uniqueness of the weak-channel crossing.
        for i in (1, 2):
            correction = sum(b*g - d*t for b, d, g, t in
                             zip(running["B"][i], running["D"][i], gauge_squares, traces))
            assert running["BETA"][i] + correction/loop < 0

    gauge["evolve"](low, logs, q(33), 8192, diagonal_quarks=True, observe=observe)
    # x_norm was obtained by adding bounds on the increments to its initial 1.
    defect = x_norm - 1
    assert 0 < defect < q(".023456")
    return defect


def velocity_sum(index, time):
    common = gauge["target"](index, time)
    terminal = [common*r for r in ratios]
    traces = [arb(0).union(arb(1))]*3
    total = arb(0)
    for i, r in enumerate(ratios):
        correction = sum(b/x - d*t for b, d, x, t in
                         zip(running["B"][i], running["D"][i], terminal, traces))
        velocity = 2*common - 2*(running["BETA"][i] + correction/loop)/(loop*r)
        assert velocity > 0
        total += velocity
    return total


def verify():
    ctx.prec = 144
    gauge_factors = q(".999").union(q("1.001"))
    yukawa_factors = q(".99").union(q("1.01"))
    low = [1/(c*gauge_factors)**2 for c in gauge["centers"]]
    logs = [2*(c*yukawa_factors).log() for c in gauge["yukawa_centers"]]
    defect = gauge_derivative_defect(low, logs)

    first, last = q("26.81235"), q("26.81262")
    inner = first.union(last)
    target = gauge["target"](11488, inner)
    margins = []
    for i, r in enumerate(ratios):
        for factor, sign in ((q(".999"), 1), (q("1.001"), -1)):
            face = low.copy()
            face[i] = 1/(gauge["centers"][i]*factor)**2
            end, _, _ = gauge["evolve"](face, logs, inner, 8192, diagonal_quarks=True)
            margins.append((sign*(end[i]/r - target)).lower())
    assert min(margins) > q(".0001488166")

    left, right = q("26.859461"), q("26.859504")
    for time, i, sign in ((left, 1, -1), (right, 2, 1)):
        end, _, _ = gauge["evolve"](low, logs, time, 8192, diagonal_quarks=True)
        assert sign*(end[i]/ratios[i] - gauge["target"](10875, time)) > 0

    ratio = ((last-first)/(right-left) * velocity_sum(11488, inner)
             / velocity_sum(10875, left.union(right)) * ((1-defect)/(1+defect))**3)
    w_other, w_top = [spectral["weight"](acb.zeta_zero(n).imag) for n in (11488, 10875)]
    weighted = ratio*w_other/w_top
    assert 0 < w_other < w_top
    assert ratio > q("5.44765") and weighted > q("4.67295")
    print("Gauge Jacobian row-norm defect <", defect.upper().str(9),
          "minimum supporting face margin", min(margins).str(9))
    print("11488/10875 full low-input coefficient >", ratio.lower().str(9),
          "spectrally weighted ratio >", weighted.lower().str(9))


if __name__ == "__main__":
    verify()
