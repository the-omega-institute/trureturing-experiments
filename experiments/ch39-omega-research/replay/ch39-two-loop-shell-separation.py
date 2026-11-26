"""Interval enclosure of a two-loop SM gauge trajectory against the local Weyl face bound."""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ctx.prec = 160
ROOT = Path(__file__).resolve().parent
local = run_path(str(ROOT/'ch39-local-weyl-coordination.py'))


def q(text):
    value = F(text)
    return arb(value.numerator)/value.denominator


def matrix(entries):
    return [list(map(q, row.split())) for row in entries]


B = matrix(['199/50 27/10 44/5', '9/10 35/6 12', '11/10 9/2 -26'])
D = matrix(['17/10 1/2 3/2', '3/2 3/2 1/2', '2 2 0'])
C = matrix(['17/20 9/4 8', '1/4 9/4 8', '9/4 9/4 0'])
BETA = list(map(q, ['41/10', '-19/6', '-7']))
LOOP = 16*arb.pi()**2


def hull_inputs(entries, factor):
    return [q(x)*q(1-factor).union(q(1+factor)) for x in entries.split()]


def integrate(gauge, steps, observe=None, *, subdivisions=64):
    start = [1/g**2 for g in gauge]
    yukawa = hull_inputs('.00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101', F(1, 100))
    y = [2*v.log() for v in yukawa]
    correction = [arb(0) for _ in range(3)]
    result = None
    step = q(F(1, subdivisions))
    for j in range(steps):
        time = j*step+arb(0).union(step)
        wide_x = [u+arb(0, '0.0002') for u in correction]
        wide_y = [u+arb(0, '0.004') for u in y]
        residual = [a+b for a, b in zip(start, wide_x)]
        inverse = [a-2*b*time/LOOP for a, b in zip(residual, BETA)]
        assert all(g > q('1/2') for g in inverse)
        if j >= 23*subdivisions:
            assert 4*arb.pi()*inverse[1] > 42
        gs = [1/v for v in inverse]
        ys = [v.exp() for v in wide_y]
        tr = [sum(ys[i:i+3]) for i in (0, 3, 6)]
        total = 3*(tr[0]+tr[1])+tr[2]
        dx = [-2*sum(b*g-d*h for b, g, d, h in zip(br, gs, dr, tr))/LOOP**2 for br, dr in zip(B, D)]
        if observe is not None:
            observe(j, inverse, [v-2*b/LOOP for v, b in zip(dx, BETA)])
        dy = []
        for f, row in enumerate(C):
            other = arb(0).union(tr[1-f]) if f < 2 else arb(0)
            gy = sum(c*g for c, g in zip(row, gs))
            dy.extend((3*(v-other)+2*(total-gy))/LOOP for v in ys[3*f:3*f+3])
        assert all(abs(step*d) < q('.0002') for d in dx)
        assert all(abs(step*d) < q('.004') for d in dy)
        # G x beta = (G + 2 beta*time/LOOP) x beta. Remove the entire one-loop term symbolically.
        axis0 = residual[1]*BETA[2]-residual[2]*BETA[1]
        axis1 = residual[2]*BETA[0]-residual[0]*BETA[2]
        assert axis0 < 0
        coordinate = -axis1/axis0
        result = coordinate if result is None else result.union(coordinate)
        correction = [v+step*d for v, d in zip(correction, dx)]
        y = [v+step*d for v, d in zip(y, dy)]
    end = [a+r-2*b*(steps*step)/LOOP for a, r, b in zip(start, correction, BETA)]
    return result, end


def enclosure(observe=None):
    gauge = hull_inputs('.467453 .63811 1.0583', F(1, 1000))
    return integrate(gauge, 33*64, observe)[0]


def verify():
    coordinate = enclosure()
    lo, hi = F('2.4727'), F('2.5416')
    assert q(lo) < coordinate < q(hi)
    normals = []
    for u in (lo, hi):
        n = (F(-1), u, (-F(41, 10)-F(19, 6)*u)/7)
        den = lcm(*(x.denominator for x in n))
        normals.append(tuple(int(den*x) for x in n))
    states = {local['TRIVIAL']: 1}
    tables = {}
    for d in range(1, 18):
        states = local['multiply_ledger'](states)
        if d in (7, 13, 17):
            tables[d] = states
    for degrees in ((17, 17, 17), (7, 13, 13)):
        gap, _, critical = local['bounds'](tables, degrees, normals)
        print(degrees, float(gap), float(critical/(4-30*critical)))
        assert gap > (F('0.000319') if degrees[0] == 17 else F('0.03070'))
        weak_gap, _, _ = local['bounds'](tables, degrees, [(0, -1, 0)])
        upper_weak = 2-weak_gap
        threshold = 42*gap/(upper_weak+gap)
        assert threshold > (F('0.0129') if degrees[0] == 17 else F('1.24'))
        print('Minimum ordinary threshold in inverse-alpha units:', float(threshold))
    print('Two-loop normal coordinate in [2.4727, 2.5416].')


if __name__ == '__main__':
    verify()
