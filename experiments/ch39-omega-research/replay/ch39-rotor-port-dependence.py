"""Certify port dependence of the three- and four-step signed return potential."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def ball(value):
    return arb(value.numerator) / value.denominator


def mean_coefficients(rows):
    return tuple(sum(column, Fraction()) / len(rows) for column in zip(*rows))


def verify():
    ctx.prec = 192
    shell = run_path(str(Path(__file__).with_name("ch39-rotor-shell-returns.py")))
    nodes, adj, root = shell["charged_shell"]()
    coefficients = []
    for index in range(len(nodes)):
        _, _, three, four = shell["returns"](adj, index)
        assert all(abs(k) <= 1 for k in three) and all(abs(k) <= 2 for k in four)
        assert all(three[-k] == v for k, v in three.items())
        assert all(four[-k] == v for k, v in four.items())
        coefficients.append((three.get(1, Fraction()), four.get(1, Fraction()),
                             four.get(2, Fraction())))

    def cosine(zero):
        return 1 - 1 / (2 * (zero.imag**2 + arb(1) / 4))

    first = cosine(acb.zeta_zero(1))
    x = [cosine(z) for z in acb.zeta_zeros(10874, 3)]
    assert first < x[0] < x[1] < x[2] < 1
    parameter = Fraction(7460124, 259581545)

    def half_slope(row, position):
        b, a, c = row
        return ball(a - parameter*b) + 4*ball(c)*position

    midpoints = [(u + v)/2 for u, v in zip(x, x[1:])]
    counts = Counter()
    intervals = []
    for row, multiplicity in Counter(coefficients[:19]).items():
        b, a, c = row
        if row == (0, 0, 0):
            counts["flat"] += multiplicity
            continue
        assert b > 0 and c > 0
        intervals.append(tuple((ball(a) + 4*ball(c)*m)/ball(b) for m in midpoints))
        if half_slope(row, first) > 0:
            counts["first"] += multiplicity
        elif half_slope(row, 1) < 0:
            counts["decreasing"] += multiplicity
        else:
            assert row == coefficients[root]
            assert half_slope(row, midpoints[0]) < 0 < half_slope(row, midpoints[1])
            counts["10875"] += multiplicity
    assert counts == {"flat": 7, "first": 6, "decreasing": 3, "10875": 3}
    intervals.sort(key=lambda pair: float(pair[0]))
    assert all(left[1] < right[0] for left, right in zip(intervals, intervals[1:]))
    assert all(c == 0 for _, _, c in coefficients[19:])

    old = mean_coefficients(coefficients[:19])
    full = mean_coefficients(coefficients)
    assert half_slope(old, 1) < 0 and half_slope(full, first) > 0

    for name, average, side, lower, upper in [
        ("19-point mean", old, 1, "3.37115e-15", "3.37116e-15"),
        ("76-point mean", full, 0, "2.02103e-15", "2.02104e-15"),
    ]:
        local = [half_slope(coefficients[root], m) for m in midpoints]
        mixed = [half_slope(average, m) for m in midpoints]
        assert local[0] < 0 < local[1]
        assert local[1-side] * mixed[1-side] > 0
        assert local[side] * mixed[side] < 0
        crossing = local[side] / (local[side] - mixed[side])
        assert arb(lower) < crossing < arb(upper)
        print(name, "first tie at mixture fraction", crossing.str(17))
    print("Original ports by minimum:", dict(counts))


if __name__ == "__main__":
    verify()
