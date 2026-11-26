"""Two finite conjugate color/weak rotations separate zero labels despite a central disk factor."""
from fractions import Fraction as F
from json import loads
from pathlib import Path

ROOT = Path(__file__).resolve().parent
T = F(1, 15)
U = 1+T
GRID = 1 << 128
DEGREE = 28


def down(x):
    return (x.numerator*GRID)//x.denominator


def up(x):
    return -down(-x)


COEFFICIENTS = [(down((U**k-1)/k), up((U**k-1)/k)) for k in range(1, DEGREE+1)]


def series(x):
    lo_x, hi_x = down(x), up(x)
    lo = hi = 0
    for a, b in reversed(COEFFICIENTS):
        lo = a+lo*lo_x//GRID
        hi = b-((-hi*hi_x)//GRID)
    # For x<=1/7, one finite rational bound works for all 3200 rotations.
    assert 0 < x <= F(1, 7)
    return lo, hi+TAIL


TAIL = up(U*(U/7)**DEGREE/((DEGREE+1)*(1-U/7)))


def signature(x):
    a, b = series(x)
    c, d = series(2*x)
    return F(2*c, b), F(2*d, a)


def verify():
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    rows = []
    for cell, *_ in prefix:
        lo, hi = (F(k, 1 << 400) for k in (cell, cell+1))
        offset = F(1, 4)+30+2*T
        x_min, x_max = 16/(hi**2+offset), 16/(lo**2+offset)
        rows.append((signature(x_min)[0], signature(x_max)[1]))
    assert len(rows) == 3200 and rows[0][1] < F('2.2')
    gaps = [rows[i][0]-rows[i+1][1] for i in range(len(rows)-1)]
    worst = min(range(len(gaps)), key=gaps.__getitem__)
    assert worst+1 == 3150 and min(gaps) > F('2.0262e-10')
    # For N >= the microscopic count, bound the effect of relative error in all four weights.
    error = F('1e-10')
    log_error = 2*error/(1-error)
    assert 144196081*T*x_min > 11
    signature_error = (1+rows[0][1])*log_error/(11-log_error)
    assert 2*signature_error < min(gaps)
    for bits in range(64):
        scale = 1 << bits
        codes = [((lo-2)*scale//1, (hi-2)*scale//1) for lo, hi in rows]
        assert all(a == b for a, b in codes)
        if len(set(codes)) == len(codes):
            coarse = {}
            for j, (code, _) in enumerate(codes, 1):
                coarse.setdefault(code//2, []).append(j)
            collision = next(v for v in coarse.values() if len(v) > 1)
            break
    else:
        raise AssertionError('Insufficient binary precision')
    assert bits == 32 and collision == [2915, 2916]
    print('Bits for the downward-rounded signature, previous collision:', bits, collision)
    print('Smallest separation, labels, last-to-first signature endpoints:', float(min(gaps)),
          (worst+1, worst+2), float(rows[-1][0]), float(rows[0][1]))


if __name__ == '__main__':
    verify()
