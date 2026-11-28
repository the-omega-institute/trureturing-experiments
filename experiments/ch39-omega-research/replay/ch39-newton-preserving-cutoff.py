"""Certified log-cutoff radii for a vacuum sign change on the level K=1."""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


native = run_path(str(Path(__file__).with_name("ch39-cutoff-scale-stationarity.py")))
# Zero-based weighted-median pivot and two movable nodes (driver, compensator).
CASES = {
    (11661, 1): (2, 3, 0),
    (11934, 1): (2, 3, 1),
    (12207, 1): (2, 3, 0),
    (11661, 2): (3, 4, 2),
}


def slopes(rows, xs, weights, radius):
    scale = (-radius).exp().union(radius.exp())
    a, b = [], []
    for x, weight in zip(xs, weights, strict=True):
        band = x*scale
        assert band > 0 and band < 1
        _, _, phi, psi = native["heat"](rows, band, 1)
        a.append(weight*phi[0]/(2*arb.pi()))
        b.append(-weight*psi[0]/(4*arb.pi()))
    return a, b


def descent(coefficient):
    assert not coefficient.contains(0)
    return -1 if coefficient > 0 else 1


def compensated_drop(a, b, pivot, velocities):
    correction = -sum(value*speed for value, speed in zip(a, velocities, strict=True))/a[pivot]
    assert abs(correction) < 1
    return -sum(value*speed for value, speed in zip(b, velocities, strict=True))-b[pivot]*correction


def verify():
    spectra = native["spectra"]()
    ctx.prec = 512
    c = arb("1e-122")
    for (number, order), (pivot, driver, compensator) in CASES.items():
        rows = spectra[number]
        table = native["FIRST"] if order == 1 else native["SECOND"]
        fixed, lo, hi = table[number]
        xs, weights, _, _ = native["construct"](rows, fixed, lo, hi, order)
        a0, b0 = slopes(rows, xs, weights, arb(0))
        multiplier = b0[pivot]/a0[pivot]
        movements = [descent(b-multiplier*a) if j != pivot else 0
                     for j, (a, b) in enumerate(zip(a0, b0, strict=True))]
        tangent = compensated_drop(a0, b0, pivot, movements)
        # Reuse a uniform box strictly larger than every certified flip.
        radius = arb("1e-113")
        a, b = slopes(rows, xs, weights, radius)
        drop = compensated_drop(a, b, pivot, movements)
        upper_slope = sum(abs(y-multiplier*x) for x, y in zip(a, b, strict=True))
        lower_radius = (c/upper_slope).lower()
        upper_radius = (arb("1.0001")*c/drop).upper()
        assert drop > 0 and 0 < lower_radius and upper_radius < radius
        assert upper_radius/lower_radius < arb("1.001")
        pair_speed = [0]*len(xs)
        pair_speed[driver] = descent(b0[driver]-b0[compensator]*a0[driver]/a0[compensator])
        pair_drop = compensated_drop(a, b, compensator, pair_speed)
        assert pair_drop > 0
        pair_radius = (arb("1.0001")*c/pair_drop).upper()
        assert pair_radius < radius
        print(number, order, "sensitivity", tangent.str(7), "radius", lower_radius.str(8),
              upper_radius.str(8), "pair upper", pair_radius.str(8))


if __name__ == "__main__":
    verify()
