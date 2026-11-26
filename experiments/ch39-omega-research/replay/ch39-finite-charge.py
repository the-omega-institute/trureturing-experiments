"""Subbands of the charged shell and finite right spectra matching its first two heat moments."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / 'ch39-slope-matched-one-loop.py'))
q = model['q']
h = model['h']
d = model['d']
pi = arb.pi()
scale = q(400)**2
cap = (-q(1)).exp()
low_y = [energy / scale for energy in model['xs']]
low_sum = sum((-y).exp() for y in low_y)
low_moment = sum(y * (-y).exp() for y in low_y)

def poles_below(y):
    # Only the first harmonic can enter the two test intervals.
    energy = scale * y
    assert energy < 4*d
    count = 0
    max_row_term = q(0)
    for mult, factor in model['factors']:
        k = pi * (energy / factor).sqrt()
        max_row_term = max(max_row_term, abs(k / k.sin()).upper())
        if factor < energy:
            count += mult
        else:
            assert factor > energy
    assert 959990 - energy - 36*max_row_term > 0
    return count

def band(ylo, yhi, last_pole, next_cut):
    # The lower edge uses the continuum count followed by the lumped-mass comparison.
    low, high, pole, gap = map(q, (ylo, yhi, last_pole, next_cut))
    beta = model['beta']
    def continuum_test(finite_y):
        return finite_y * (1+10*beta) / (1-beta*scale*finite_y)
    first = 19 + poles_below(continuum_test(low))
    last = 18 + poles_below(continuum_test(gap))
    # Trial space: static extensions of 18 low endpoint components plus Dirichlet modes.
    prior = [(n, factor) for n, factor in model['factors'] if factor/scale < low]
    assert first == 19+sum(n for n, _ in prior)
    included = sum(n for n, factor in model['factors'] if factor/scale < pole)
    assert last == 18+included
    static_mass = 9*pi**2/d
    dirichlet_bound = scale*pole/(1-h*h)**2
    previous = max((factor for _, factor in prior), default=q(0))/(1-h*h)**2
    test_low = scale*low
    assert test_low > previous and (test_low-30)*(test_low-previous) > test_low**2*static_mass
    test = scale*high
    assert test > 30 and test > dirichlet_bound
    assert (test-30)*(test-dirichlet_bound) > test**2*static_mass
    print('ordered band', first, last, 'mass squared / Q squared', low.str(8), high.str(8))
    return low.union(high), first, last

first_band, start, end = band('1.2710', '1.273', '1.27165', '1.31')
assert (start, end) == (19, 30)
selected_band, start, end = band('2.150', '2.158', '2.1550', '2.162')
assert end-start+1 == 12

def entropy(w):
    return -w*w.log()

def heat_gap(number, y, s, high_right):
    right_y = (-2*s).exp()
    right_weight = 18*(-right_y).exp()
    w = low_sum - right_weight + number*(-y).exp()
    moment = low_moment - right_y*right_weight + number*y*(-y).exp()
    packed = int((w/cap).lower().floor().unique_fmpz())
    assert packed == int((w/cap).upper().floor().unique_fmpz())
    remainder = w-packed*cap
    assert remainder > 0 and remainder < cap
    minimum = packed*cap + entropy(remainder)
    uniform = w*(q(high_right)/w).log()
    return moment-minimum, uniform-moment, w, moment

def check(number, y, high_right, feasible):
    margins, room = [], []
    # Outward rounding covers unequal eigenvalues anywhere inside the same band.
    for j in range(55):
        s = (q(2550+j)/1000).union(q(2551+j)/1000)
        gap, upper_gap, w, _ = heat_gap(number, y, s, high_right)
        if feasible:
            assert w < high_right*cap and upper_gap > 0 and gap > 0
        else:
            assert gap < 0
        margins.append(gap)
        room.append(upper_gap)
    print('high left/right',number,high_right,'moment minus packed entropy',
          min(g.lower() for g in margins).str(9),max(g.upper() for g in margins).str(9),
          'uniform entropy minus moment lower',min(g.lower() for g in room).str(9))

for n in range(1,6):
    check(n, first_band, n, n == 5)
check(12, first_band, 12, True)
check(2, selected_band, 2, True)
right_y = (-2*q('2.55').union(q('2.605'))).exp()
assert low_sum-18*(-right_y).exp() > 0
assert low_moment-18*right_y*(-right_y).exp() < 0
