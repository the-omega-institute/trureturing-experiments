"""Real one-loop separation from a local bound on the three stars of an observed triangle."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import gcd, lcm, prod
from pathlib import Path
import json

TRIVIAL = (0, 0, 0, 0)
UNIT = 200
GRID = 1 << 100
BARE = F(4, UNIT+30)


def multiply_ledger(states):
    out = defaultdict(int)
    for (y, n, p, q), coeff in states.items():
        color = [(p+1, q)]
        anti = [(p, q+1)]
        if p:
            color.append((p-1, q+1))
            anti.append((p-1, q))
        if q:
            color.append((p, q-1))
            anti.append((p+1, q-1))
        out[y, n, p, q] += UNIT*coeff
        for charge in (y+6, y-6):
            out[charge, n, p, q] += coeff
        for weak in ([n+1, n-1] if n else [1]):
            for charge in (y+3, y-3):
                out[charge, weak, p, q] += coeff
            for a, b in color:
                out[y+1, weak, a, b] += coeff
            for a, b in anti:
                out[y-1, weak, a, b] += coeff
        for a, b in color:
            out[y+4, n, a, b] += coeff
            out[y-2, n, a, b] += coeff
        for a, b in anti:
            out[y-4, n, a, b] += coeff
            out[y+2, n, a, b] += coeff
    return out


def normal(em, sw, strong):
    start = (3*em*(1-sw)/5, em*sw, 1/strong)
    beta = (F(41, 10), F(-19, 6), F(-7))
    v = [start[(j+1)%3]*beta[(j+2)%3]-start[(j+2)%3]*beta[(j+1)%3] for j in range(3)]
    den = lcm(*(q.denominator for q in v))
    integers = [int(q*den) for q in v]
    return tuple(n//gcd(*integers) for n in integers)


def bounds(tables, degrees, normals):
    data = [tables[h] for h in degrees]
    scale = UNIT**sum(degrees)
    vacuum = prod(2*UNIT**h-table[TRIVIAL] for h, table in zip(degrees, data))
    assert all(states[TRIVIAL] < 2*UNIT**h for h, states in zip(degrees, data))
    low = vacuum*GRID//scale
    err = [0]*len(normals)
    for (y, n, p, q), first in data[0].items():
        if (y, n, p, q) == TRIVIAL:
            continue
        dim = (n+1)*(p+1)*(q+1)*(p+q+2)//2
        key = (y, n, p, q)
        numerator = first*data[1].get(key, 0)*data[2].get(key, 0)
        value = -((-numerator*GRID)//(scale*dim))
        low -= value
        # 120 times the mean square of each unified generator on this representation.
        charges = (2*y*y, 10*n*(n+2), 5*(p*p+p*q+q*q+3*p+3*q))
        for i, (a, b, c) in enumerate(normals):
            err[i] += value*abs(a*charges[0]+b*charges[1]+c*charges[2])
    assert low > 0
    gaps = [(-sum(axis)-F(e, 120*low)/BARE)/sum(map(abs, axis))
            for e, axis in zip(err, normals)]
    critical = max(F(e, 120*low)/(-sum(axis)) for e, axis in zip(err, normals))
    return min(gaps), F(low, GRID), critical


def verify():
    prefix = json.loads(Path(__file__).with_name('ch39-tick-separators.json').read_text())['prefix']
    cell = prefix[0][0]
    assert F(1, 4)+F(cell, 1 << 400)**2 > UNIT
    normals = [normal(a, s, g) for a, s, g in product(
        (F('127.90'), F('127.96')), (F('.23125'), F('.23133')), (F('.117'), F('.119')))]
    states = {TRIVIAL: 1}
    tables = {}
    for degree in range(1, 18):
        states = multiply_ledger(states)
        if degree in (7, 13, 17):
            tables[degree] = states
            assert sum(v*(n+1)*(p+1)*(q+1)*(p+q+2)//2 for (y, n, p, q), v in states.items()) == (UNIT+30)**degree
    for degrees in ((17, 17, 17), (7, 13, 13)):
        gap, vacuum, critical = bounds(tables, degrees, normals)
        cutoff = critical/(4-30*critical)
        count = 1 if degrees[0] == 17 else 6
        last = F(prefix[count-1][0]+1, 1 << 400)
        following = F(prefix[count][0], 1 << 400)
        assert 1/(F(1, 4)+last**2) > cutoff > 1/(F(1, 4)+following**2)
        assert gap > (F('.00622') if degrees[0] == 17 else F('.03'))
        print('Exterior face counts', degrees, 'denominator lower bound', float(vacuum),
              'separation in bare stiffness units', float(gap), 'cutoff for arbitrary neighbors', float(cutoff))
    h = {k: 2**k+2**(5-k)-5 for k in range(1, 5)}
    assert all(all(a <= b for a, b in zip(sorted((h[i], h[j], h[i+j])), (7, 13, 13)))
               for i in h for j in h if i+j <= 4)


if __name__ == '__main__':
    verify()
