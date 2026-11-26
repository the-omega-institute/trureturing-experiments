"""Full SU(5) Haar edges, with a rank-two frame inserted in each based face trace."""
from fractions import Fraction as F
from json import loads
from math import prod
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
model = run_path(str(ROOT/'ch39-uniform-weak-center.py'))
q = model['q']
t = model['AMPLITUDE']


def exterior(strength):
    # One insertion for each of 5, 5-bar, 10, 10-bar.
    first = [h*(1+t)*strength for h in (7, 13, 13)]
    x = [30*v for v in first]
    assert max(x) < 1
    masses = [v/(1-v) for v in x]
    remainders = [v*v/(2*(1-v/3)) for v in x]
    moments = [2*v*v/(1-3*v/4) for v in x]

    def rest(values):
        return sum(values[i]*masses[(i+1)%3]*masses[(i+2)%3] for i in range(3))/5**4

    dims = (5, 5, 10, 10)
    indices = (F(1, 2), F(1, 2), F(3, 2), F(3, 2))
    leading_mass = prod(first)*sum(F(1, d) for d in dims)
    leading_moment = prod(first)*sum(b/d**2 for b, d in zip(indices, dims))
    vacuum = prod(1-v for v in remainders)-leading_mass-rest(remainders)
    assert 0 < vacuum < 1
    return (leading_moment+rest(moments))/vacuum


def verify():
    cells = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    color = 4*model['arb'].pi()*q(model['COLOR'])
    unit = 1 << 40
    limits = {46: '.000243', 320: '.0000000496', 1000: '.000000000423', 3200: '.00000000000588'}
    for j, limit in limits.items():
        cell = cells[j-1][0]
        lo, hi = (F(k, 1 << 400) for k in (cell, cell+1))
        min_s, max_s = 1/(hi**2+F(1, 4)), 1/(lo**2+F(1, 4))
        bare_lower = 4*min_s/(1+(30+2*t)*min_s)
        relative = exterior(max_s)/bare_lower
        assert 4*(1+t+relative) < 5
        desired = color*(q(lo**2+F(1, 4)).union(q(hi**2+F(1, 4)))+q(30+2*t))/4
        floor = (unit*desired).floor().unique_fmpz()
        assert floor is not None
        bound = color*q(relative)+q(5*max_s/unit)
        assert bound < q(limit)
        if j == 3200:
            assert int(floor) == 158545267438980033232
        print('Full SU(5) direction bound at zero:', j, bound.str(10))


if __name__ == '__main__':
    verify()
