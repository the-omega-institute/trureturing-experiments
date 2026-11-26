"""Integer certificate for reconstruction from 12 offsets of four oriented phase channels."""
from fractions import Fraction as Q
import json
from pathlib import Path

BITS = 144
OFFSETS = (0,4,7,11,15,19,22,26,30,34,37,41)


def turn(chart, den):
    sign, n = chart
    return Q(sign*(den**2 - n*n), den**2 + n*n), Q(2*sign*n*den, den**2 + n*n)


def matrix(table):
    den = 2**table["bits"]
    rows = []
    for j in OFFSETS:
        block = [[] for _ in range(4)]
        for _, h, w, a, b in table["rows"]:
            ca, sa = turn(a, den)
            cb, sb = turn(b, den)
            sines = (sa, sb, ca*sb-sa*cb, ca*sb+sa*cb)
            mass = Q(w*h*h, den*(h*h + (j*den)**2))
            for row, s in zip(block, sines):
                row.append(mass*s*Q(7,68))
        rows.extend(block)
    return rows[:46]


def main():
    root = Path(__file__).resolve().parent
    original = matrix(json.loads((root/"ch39-five-tick-record.json").read_text()))
    inv = json.loads((root/"ch39-oriented-height-inverse.json").read_text())
    den = 2**BITS
    assert len(inv) == len(original) == 46 and all(len(row) == 46 for row in inv+original)
    grid = [[x*den//1 for x in row] for row in original]
    # grid/den encloses each exact entry to within 1/den.
    # The true response from the 288-bit source cells also fits within an
    # additional 1/den per entry (verified by the analytic origin script).
    max_error = 0
    max_inverse = 0
    for i, row in enumerate(inv):
        size = sum(map(abs, row))
        err = 0
        for j in range(46):
            residual = sum(row[k]*grid[k][j] for k in range(46)) - int(i == j)*den**2
            err += abs(residual) + 2*size
        max_error = max(max_error, err)
        max_inverse = max(max_inverse, size)
    assert 2*max_error < den**2
    assert max_inverse*den < 2**73*(den**2-max_error)
    print("At the twelve offsets, input errors <= 2^-96 give coefficient errors < 2^-23.")
    print("The integer inequality covers any matrix within 2^-144 per entry of the dyadic-clock matrix.")


if __name__ == "__main__":
    main()

