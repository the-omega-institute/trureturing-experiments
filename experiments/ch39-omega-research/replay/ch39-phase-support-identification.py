"""Identify a 50-site positive spectrum with 68 time rows and one shell observable."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
time = run_path(str(root / "ch39-time-support-identification.py"))
base, B, rounded, multiply = [time[k] for k in ("base", "B", "rounded", "multiply")]
controls = [33, 35, 46, 49, 59, 63, 82, 85, 87, 113, 121, 128, 196, 211,
            437, 456, 474, 2314]
coefficients = tuple(Fraction(n, 10**8) for n in (-64408846, -68137085, -13409845, 100000000))
# Nonnegative slack at the interpolation sites, in normalized residue units.
slack = dict(zip(controls, (0, 2160037, 642866, 1280588, 697409, 1492007, 105392, 0,
                            2953752, 567498, 1516891, 1836234, 565740, 889832,
                            6299230072, 4242790287, 3206682140, 461918905074)))
offsets = [*base["OFFSETS"], *range(43, 67)]
unit = 1 << 185


def integer_cell(value):
    return value.lo * unit // 1, -((-value.hi * unit) // 1)


def height_cell(floor):
    return rounded(B(Fraction(1 << 800, (floor + 1)**2), Fraction(1 << 800, floor**2)))


def shell_values(ids):
    matrix, _ = base["observations"](ids)
    rows = matrix[-4:]
    return [sum((B(c) * B(Fraction(lo, base["SCALE"]), Fraction(hi, base["SCALE"]))
                 for c, (lo, hi) in zip(coefficients, column)), B(0))
            for column in zip(*rows)]


def interpolation_data(ids, heights):
    nodes = [heights[i] for i in ids]
    betas = []
    for k, x in enumerate(nodes):
        numerator = multiply(*(B(1) + B(j*j)*x for j in offsets if j))
        denominator = multiply(*(x-y for m, y in enumerate(nodes) if k != m))
        betas.append(rounded(numerator / denominator))
    scale = B(max(max(abs(v.lo), abs(v.hi)) for v in betas))
    data = []
    for i, x, b, f in zip(ids, nodes, betas, shell_values(ids)):
        assert b.lo > 0 or b.hi < 0
        charge = Fraction(slack.get(i, 0), 10**12) * (1 if b.lo > 0 else -1)
        data.append((integer_cell(x), integer_cell(b/scale), integer_cell(f), integer_cell(B(charge))))
    return data


def residual_floor(cell, value, data):
    # Multiply the interpolation difference by W/|H| and a positive scale.
    sign = (-1)**sum(cell[1] < x[0] for x, *_ in data)
    lower = 0
    for (xl, xh), (bl, bh), (fl, fh), charge in data:
        dl, dh = cell[0]-xh, cell[1]-xl
        assert dl*dh > 0
        yl, yh = value-fh, value-fl
        products = (bl*yl, bl*yh, bh*yl, bh*yh)
        low, high = min(products) + unit*charge[0], max(products) + unit*charge[1]
        if sign < 0:
            low, high = -high, -low
        numerator = low if dl > 0 else high
        lower += min(numerator // dl, numerator // dh)
    return lower


def main():
    roster = set(time["roster"]()[0])
    ids = sorted(roster | set(controls))
    assert len(roster) == 50 and len(ids) == len(set(offsets)) == 68 and 0 in offsets
    assert all(v >= 0 for v in slack.values()) and not (set(slack) & roster)
    cells = base["data"]("ch39-tick-separators.json")["prefix"]
    assert len(cells) == 3200 and all(a[0]+1 < b[0] for a,b in zip(cells,cells[1:]))
    heights = {i: height_cell(row[0]) for i, row in enumerate(cells, 1)}
    data = interpolation_data(ids, heights)
    candidates = set(ids)
    lower = (-sum(map(abs, coefficients)) * unit) // 1
    exceptions = [i for i, x in heights.items()
                  if i not in candidates and residual_floor(integer_cell(x), lower, data) <= 0]
    for i, f in zip(exceptions, shell_values(exceptions)):
        assert residual_floor(integer_cell(heights[i]), integer_cell(f)[0], data) > 0
    # Cover the entire real interval below the last tabulated ordinate.
    todo = [(0, integer_cell(heights[len(cells)])[1])]
    covered = 0
    while todo:
        lo, hi = todo.pop()
        # The linear interpolation error is bounded by width^2 times the sum
        # of |numerator|/distance^3, since the second derivative has factor 2.
        loss = 0
        for (xl, _), beta, f, charge in data:
            distance = xl - hi
            assert distance > 0
            bound = max(map(abs, beta))*max(abs(lower-y) for y in f) + unit*max(map(abs, charge))
            loss += -((-((hi-lo)**2 * bound)) // distance**3)
        if min(residual_floor((t, t), lower, data) for t in (lo, hi)) > loss:
            covered += 1
        else:
            mid = (lo+hi)//2
            assert lo < mid < hi
            todo.extend([(lo, mid), (mid, hi)])
    print("68 time rows and one linear combination of the four shell means suffice.")
    print("Exact shell words needed at", len(exceptions), "additional prefix sites.")
    print("Closed intervals covering all higher zeros, even with unrestricted phases:", covered)


if __name__ == "__main__":
    main()
