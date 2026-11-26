"""Finite positive Weyl subdivisions and exponential plaquette weights at fixed total strength."""
from fractions import Fraction as F
from math import factorial, lcm, prod
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
local = run_path(str(ROOT/'ch39-local-weyl-coordination.py'))
TRIVIAL = local['TRIVIAL']
GRID = 1 << 100
DEGREE = 20


def tables(degrees, strength):
    states = {TRIVIAL: 1}
    num, unit = strength.as_integer_ratio()
    den = factorial(DEGREE)*unit**DEGREE
    out = {h: {} for h in degrees}
    for k in range(DEGREE+1):
        weights = {h: (h*num)**k*(factorial(DEGREE)//factorial(k))*unit**(DEGREE-k) for h in degrees}
        assert sum(c*(n+1)*(p+1)*(q+1)*(p+q+2)//2 for (y, n, p, q), c in states.items()) == 30**k
        for h in degrees:
            table, weight = out[h], weights[h]
            for irrep, coeff in states.items():
                table[irrep] = table.get(irrep, 0)+weight*coeff
        if k < DEGREE:
            expanded = local['multiply_ledger'](states)
            for irrep, value in states.items():
                expanded[irrep] -= 200*value
            states = {k: v for k, v in expanded.items() if v}
    return out, den


def tails(h, strength):
    x = 30*h*strength
    k = DEGREE+1
    # For the squared-generator tail, the next-to-current term ratio is x*(k+1)/k**2.
    ratio = x*(k+1)/k**2
    first = x**k/factorial(k)
    assert ratio < 1
    mass = first/(1-x/(k+1))
    moment = k*k*first/(1-ratio)
    total = sum(x**j/factorial(j) for j in range(k))+mass
    return mass, moment, total


def exterior(rows, den, degrees, normals, distortion):
    strength = (1+distortion)/F(200)
    bare_lower = F(4, 200)/(1+30*strength)
    data = [rows[h] for h in degrees]
    rest = [tails(h, strength) for h in degrees]
    low = prod(2-F(row[TRIVIAL], den)-tail[0] for row, tail in zip(data, rest))
    assert all(2-F(row[TRIVIAL], den)-tail[0] > 0 for row, tail in zip(data, rest))
    scale = den**3
    low_grid = low*GRID//1
    err = [0]*len(normals)
    for key, first in data[0].items():
        if key == TRIVIAL:
            continue
        y, n, p, q = key
        dim = (n+1)*(p+1)*(q+1)*(p+q+2)//2
        w = -((-GRID*first*data[1].get(key, 0)*data[2].get(key, 0))//(scale*dim))
        low_grid -= w
        c1, c2, c3 = 2*y*y, 10*n*(n+2), 5*(p*p+p*q+q*q+3*p+3*q)
        for j, (a, b, c) in enumerate(normals):
            err[j] += w*abs(a*c1+b*c2+c*c3)
    denominator = F(low_grid, GRID)-sum(rest[i][0]*rest[(i+1)%3][2]*rest[(i+2)%3][2] for i in range(3))
    remainder = sum(rest[i][1]*rest[(i+1)%3][2]*rest[(i+2)%3][2] for i in range(3))
    assert denominator > 0
    errors = [(F(e, 120*GRID)+sum(map(abs, axis))*remainder)/(bare_lower*denominator)
              for axis, e in zip(normals, err)]
    return errors, denominator


# Ten times the Weyl indices of Q, u, d, L, e, in hypercharge, weak, color order.
INDICES = ((1,15,10),(8,0,5),(2,0,5),(3,5,0),(6,0,0))


def bounds(rows, den, degrees, normals, distortion):
    errors, denominator = exterior(rows, den, degrees, normals, distortion)
    gaps = []
    for axis, error in zip(normals, errors):
        norm = sum(map(abs, axis))
        defect = sum(abs(sum(a*b for a, b in zip(axis, row))) for row in INDICES)
        gaps.append((-sum(axis)-distortion*defect/20-error)/norm)
    return gaps, denominator


def verify():
    normals = []
    for u in (F('2.4727'), F('2.5416')):
        axis = (F(-1), u, (-F(41, 10)-F(19, 6)*u)/7)
        scale = lcm(*(x.denominator for x in axis))
        normals.append(tuple(int(x*scale) for x in axis))
    normals.append((0, -1, 0))
    for distortion, cases in [
        (F(0), [((7, 13, 13), '1.2290'), ((16, 16, 16), '.2414')]),
        (F('.01'), [((7, 13, 13), '.9737')]),
        (F('.03'), [((7, 13, 13), '.4681')]),
        (F('.04'), [((7, 13, 13), '.2177')]),
        (F('.05'), [((7, 13, 13), None)]),
    ]:
        needed = {h for degrees, _ in cases for h in degrees}
        rows, den = tables(needed, (1+distortion)/200)
        for degrees, claimed in cases:
            gaps, vacuum = bounds(rows, den, degrees, normals, distortion)
            gap = min(gaps[:2])
            if claimed is None:
                assert gap <= 0
                print('Inconclusive exterior bound at', distortion, degrees)
                continue
            assert gap > 0 and 2-gaps[2] > 0
            displacement = 42*gap/(2-gaps[2]+gap)
            assert displacement > F(claimed)
            print('Distortion, face counts, denominator, gap, inverse-alpha distance:',
                  distortion, degrees, float(vacuum), float(gap), float(displacement))


if __name__ == '__main__':
    verify()
