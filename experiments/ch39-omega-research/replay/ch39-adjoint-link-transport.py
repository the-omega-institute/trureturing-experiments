"""Haar response of a rank-two frame transported on each SU(5) link."""
from fractions import Fraction as F
from itertools import combinations
from json import loads
from math import prod
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
model = run_path(str(ROOT/'ch39-uniform-weak-center.py'))
q = model['q']
KAPPA = F(1, 8)
C = F(158545267438980033232, 1 << 40)
ETA = KAPPA**3/24**2
BROKEN_MOMENT = F(5, 12)*ETA


def bounds(s):
    x = [h*30*(1+model['AMPLITUDE'])*s for h in (7, 13, 13)]
    weight = 1+24*KAPPA
    mass = [weight*v/(1-v) for v in x]
    neutral = [weight*v**2/(2*(1-v/3)) for v in x]
    moment = [weight*4*v/(1-F(9, 8)*v) for v in x]
    assert max(x) < F(8, 9)
    all_active = prod(mass)/5**4
    all_active_hessian = min(moment[i]*prod(mass[k] for k in range(3) if k != i) for i in range(3))/5**4

    def mixed_mass(pure):
        return KAPPA**len(pure)*prod(neutral[i]/24 for i in range(3) if i not in pure)/24

    mixed = sum(mixed_mass(pure) for size in (1, 2) for pure in combinations(range(3), size))
    # Sorting the support counts may move the middle edge of the boundary convolution.
    unbroken_mixed = max(mixed_mass((i,)) for i in range(3))
    delta = prod(1+v for v in neutral)-1+all_active+mixed
    vacuum = 1+ETA-delta
    assert vacuum > 0
    adjoint_mean_square = F(5, 24)
    unbroken = (all_active_hessian+adjoint_mean_square*unbroken_mixed)/vacuum
    broken = ((all_active_hessian+adjoint_mean_square*mixed)/vacuum
              +BROKEN_MOMENT*delta/((1+ETA)*vacuum))
    return unbroken, broken


def verify():
    t = model['AMPLITUDE']
    cell = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix'][3199][0]
    high = F(cell+1, 1 << 400)
    low = F(cell, 1 << 400)
    max_s = 1/(F(1, 4)+low**2)
    fixed_color = 4*model['arb'].pi()*q(model['COLOR'])
    bare = 4/(q(low**2+F(1, 4)).union(q(high**2+F(1, 4)))+q(30+2*t))
    calibration_error = abs(q(C)*bare-fixed_color)*q(1+t)
    preserved, broken = bounds(max_s)
    shift = C*BROKEN_MOMENT/(1+ETA)
    unbroken_error = q(C*preserved)+calibration_error
    broken_error = q(C*broken)+calibration_error
    assert unbroken_error < q('.000000849') and broken_error < q('.00000502')
    assert F('203.72686284') < shift < F('203.72686286')
    print('Unbroken inverse-alpha error:', unbroken_error.str(10))
    print('Broken extra inverse-alpha, and error:', float(shift), broken_error.str(10))

    # The edge Hessian and the largest possible bare spatial stiffness on the broken block.
    link_hessian = F(5, 12)*KAPPA/(1+KAPPA)
    assert link_hessian == F(5, 108)
    bare_cap = 4*max_s/(1+(30+2*t)*max_s)
    norm_ratio = 42*bare_cap/link_hessian
    assert norm_ratio < F('.000261') and 1/norm_ratio > 3800
    print('Gaussian plaquette-to-link operator norm upper bound:', float(norm_ratio))


if __name__ == '__main__':
    verify()
