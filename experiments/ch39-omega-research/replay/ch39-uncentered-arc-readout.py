"""Finite color/weak rotations in an arbitrary unamplified SU(5) exterior."""
from fractions import Fraction as F
from json import loads
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
central = run_path(str(ROOT/'ch39-finite-center-pair-readout.py'))
T, U = central['T'], central['U']
GRID = 1 << 96
su5 = run_path(str(ROOT/'ch39-haar-correlation-length.py'))
ZERO = su5['ZERO']


def up(x):
    return F(-((-x*GRID)//1), GRID)


def down(x):
    return F((x*GRID)//1, GRID)


def exterior(s):
    s = up(s)
    first = [h*U*s for h in (7, 13, 13)]
    mass = [up((1+30*U*s)**h-1) for h in (7, 13, 13)]
    tail = [m-30*a for m, a in zip(mass, first)]
    nonconstant = F(3, 5)*prod(first)+sum(tail[i]*mass[(i+1)%3]*mass[(i+2)%3] for i in range(3))/5**4
    denominator = (1-(30+2*T)*s)**33
    if max(tail) < 1:
        denominator = max(denominator, prod(1-r for r in tail)-nonconstant)
    return up(2*nonconstant/denominator)


def tensor_powers():
    coeffs = [{ZERO: 1}]
    for _ in range(13):
        row = {}
        for rep, value in coeffs[-1].items():
            for child in su5['tensor'](rep):
                row[child] = row.get(child, 0)+value
        coeffs.append(row)
    return coeffs


def low_height_exterior(s, powers):
    # Precompute the majorant coefficients with finite integer multiplicities.
    u = up(U*s)
    data = []
    for h in (7, 13):
        row = {}
        for k in range(h+1):
            factor = -(int((-comb(h, k)*u**k*GRID)//1))
            for rep, value in powers[k].items():
                row[rep] = row.get(rep, 0)+factor*value
        data.append(row)
    short, long = data
    dimensions_den = prod(j-i for i, j in combinations(range(5), 2))
    mass_units = 0
    for rep, v in short.items():
        if rep == ZERO:
            continue
        dim = prod(rep[i]-rep[j]+j-i for i, j in combinations(range(5), 2))//dimensions_den
        mass_units -= (-v*long[rep]**2)//(GRID**2*dim)
    mass = F(mass_units, GRID)
    a, b = 2-F(short[ZERO], GRID), 2-F(long[ZERO], GRID)
    assert a > 0 and b > 0
    vacuum = a*b*b-mass
    assert vacuum > 0
    return up(2*mass/vacuum)


def check_roots():
    for order in permutations(range(4)):
        vertex = [0]*4
        colors = [0]
        for axis in order:
            vertex[axis] += 1
            colors.append(sum(vertex) % 5)
        assert len(set(colors)) == 5
        assert all(sum(c == 0 for c in face) <= 1 for face in combinations(colors, 3))


def verify():
    check_roots()
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    single, pair = [], []
    powers = tensor_powers()
    for index, (cell, *_) in enumerate(prefix, 1):
        lo, hi = (F(k, 1 << 400) for k in (cell, cell+1))
        strength = 1/(F(1, 4)+lo**2)
        err = exterior(strength)
        if index <= 6:
            err = min(err, low_height_exterior(strength, powers))
        offset = F(1, 4)+30+2*T
        left, right = down(16/(hi**2+offset)), up(16/(lo**2+offset))
        series_lo, series_hi = central['series'](left)[0], central['series'](right)[1]
        b_lo = down(left*series_lo/central['GRID'])
        b_hi = up(right*series_hi/central['GRID'])
        single.append((b_lo-err, b_hi+err))
        vlo, vhi = central['signature'](left)[0], central['signature'](right)[1]
        if b_lo <= err:
            pair.append(None)
        else:
            noise = up((1+vhi)*err/(b_lo-err))
            pair.append((down(vlo-noise), up(vhi+noise)))
    for label, intervals in [('single', single), ('pair', pair)]:
        gaps = [(intervals[i][0]-intervals[i+1][1]) if intervals[i] and intervals[i+1] else None
                for i in range(len(intervals)-1)]
        start = max((i+2 for i, g in enumerate(gaps) if g is None or g <= 0), default=1)
        worst = min(range(start-1, len(gaps)), key=gaps.__getitem__)
        print('Readout, certified suffix start, smallest gap, adjacent labels:', label, start,
              float(gaps[worst]), (worst+1, worst+2))
        assert worst+1 == 3150
        if label == 'single':
            assert start == 1 and min(gaps) > F('6.5361e-12')
            assert gaps[0] > F('.00234')
        else:
            assert start == 160 and min(gaps[start-1:]) > F('1.4978e-10')


if __name__ == '__main__':
    verify()
