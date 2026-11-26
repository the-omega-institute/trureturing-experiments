"""Enclose an attainable gauge trajectory and exclude the only heavier competing label."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-high-zero-planck-window.py")))
running = fixed["model"]["running"]
q = fixed["q"]
ratios = list(map(q, fixed["V"]))
centers = list(map(q, fixed["model"]["GAUGE"]))
yukawa_centers = list(map(q, ".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101".split()))
loop = running["LOOP"]


def target(index, time):
    height = acb.zeta_zero(index).imag
    factor = (fixed["geo"]["BOOK_PLANCK"] / 1000)**2 * 3 * arb(5).sqrt() / (4 * arb.pi())
    return factor * (-2 * time).exp() / (height**2 + q(fixed["OFFSET"]))


def evolve(inverse, log_y, duration, steps, *, diagonal_quarks=False, observe=None, observe_flavors=None):
    step = duration / steps
    gauge_pad, yukawa_pad = abs(step) / 7, abs(step) / 3
    b, d, c = (running[k] for k in ("B", "D", "C"))
    # Integrate the weak-minus-colour projection before taking interval bounds.
    b_diff = [w / ratios[1] - s for w, s in zip(b[1], b[2])]
    d_diff = [w / ratios[1] - s for w, s in zip(d[1], d[2])]
    one_loop = [-2 * x / loop for x in running["BETA"]]
    projected = inverse[1] / ratios[1] - inverse[2]
    for _ in range(steps):
        state = [x + arb(0, gauge_pad) for x in inverse]
        assert all(x > q("1/2") for x in state)
        gauge = [1 / x for x in state]
        squares = [(y + arb(0, yukawa_pad)).exp() for y in log_y]
        traces = [sum(squares[j:j+3]) for j in (0, 3, 6)]
        if observe is not None:
            observe(step, gauge, traces)
        if observe_flavors is not None:
            observe_flavors(step, gauge, squares, traces)
        total = 3 * (traces[0] + traces[1]) + traces[2]

        def gauge_rate(linear, weights, sources):
            return linear - 2 * sum(a*g - e*t for a, g, e, t in zip(weights, gauge, sources, traces)) / loop**2

        dx = [gauge_rate(k, row, source) for k, row, source in zip(one_loop, b, d)]
        dy = []
        for family, row in enumerate(c):
            gauge_term = sum(a*g for a, g in zip(row, gauge))
            other = arb(0).union(traces[1-family]) if family < 2 else arb(0)
            for flavor in range(3):
                mixing = squares[3*(1-family)+flavor] if diagonal_quarks and family < 2 else other
                dy.append((3*(squares[3*family+flavor]-mixing) + 2*(total-gauge_term)) / loop)
        assert all(abs(step * x) < gauge_pad for x in dx)
        assert all(abs(step * y) < yukawa_pad for y in dy)
        projected += step * gauge_rate(one_loop[1] / ratios[1] - one_loop[2], b_diff, d_diff)
        inverse = [x + step * u for x, u in zip(inverse, dx)]
        log_y = [y + step * v for y, v in zip(log_y, dy)]
    return inverse, log_y, projected


def witness():
    time = q("26.859482438")
    common = target(10875, time)
    y_factor = q(".99295").union(q(".99305"))
    logs = [2*(v*y_factor).log() for v in yukawa_centers]
    faces = [(q(".99999"), q("1.00004")), (q(".999"), q(".999005")), (q("1.000985"), q("1.001"))]
    full = [a.union(b) for a, b in faces]
    margins = []
    for i, pair in enumerate(faces):
        for side, endpoint in enumerate(pair):
            box = full.copy()
            box[i] = endpoint
            initial = [1/(v*x)**2 for v, x in zip(centers, box)]
            end, _, _ = evolve(initial, logs, time, 24576)
            residual = end[i] / ratios[i] - common
            margin = residual if side == 0 else -residual
            assert margin > q(".000002")
            margins.append(margin.lower())
    height = acb.zeta_zero(10875).imag
    curvature = 1 + (1 + (height / 312)**2).sqrt()
    assert q("35.05643642525") < curvature < q("35.05643642526")
    print("10875, minimum inward face margin", min(margins).str(9),
          "conditional curvature", curvature.str(16))


def exclusion():
    time = q(fixed["L"] + Fraction("0.04860"))
    high = target(10854, time)
    full = q(".999").union(q("1.001"))
    # Subdivide the two input coordinates that set the upper matching time.
    weak_low, weak_high = q(".999").union(q(".9992")), q(".9992").union(q("1.001"))
    color_low, color_high = q(".999").union(q("1.0008")), q("1.0008").union(q("1.001"))
    yukawa_factor = q(".99").union(q("1.01"))
    logs = [2*(v*yukawa_factor).log() for v in yukawa_centers]
    for weak, color in [(weak_high, full), (weak_low, color_low), (weak_low, color_high)]:
        initial = [1/(v*f)**2 for v, f in zip(centers, [full, weak, color])]
        end, _, difference = evolve(initial, logs, time, 8192)
        assert high - end[1] / ratios[1] > q(".00049")
        assert difference < q("-.000043")
        print("10854, weak shortfall and weak-minus-colour upper bound",
              (high-end[1]/ratios[1]).lower().str(8), difference.upper().str(8))

    def monotonic(step, inverse, deriv):
        if step >= 23*64:
            assert deriv[1] / ratios[1] > 0
            assert deriv[1] / ratios[1] < deriv[2]

    running["integrate"]([v*full for v in centers], 33*64, monotonic)


if __name__ == "__main__":
    ctx.prec = 160
    witness()
    exclusion()
