"""One finite weak-center choice and a simultaneous inverse-coupling bound for every face."""
from fractions import Fraction as F
from json import loads
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
local = run_path(str(ROOT/'ch39-exponential-face-stiffness.py'))
running = run_path(str(ROOT/'ch39-two-loop-shell-separation.py'))
q = running['q']
LOOP, BETA, B, D, C = [running[key] for key in ('LOOP', 'BETA', 'B', 'D', 'C')]
AMPLITUDE = F(1, 15)
LOG_SCALE = F(429, 16)
COLOR = F(13505, 4096)
YUKAWA_LOG_INTEGERS = (-416349, -212338, -24818, -389713, -291793, -165371, -419761, -244443, -151629)


def backward_input():
    bare = (1-F(3, 5)*AMPLITUDE, 1+AMPLITUDE, 1)
    inverse_high = [q(COLOR*x) for x in bare]
    logs = [q(F(k, 16384)) for k in YUKAWA_LOG_INTEGERS]
    correction = [arb(0) for _ in range(3)]
    step = F(1, 256)
    steps = LOG_SCALE/step
    assert steps.denominator == 1
    h = q(step)
    for j in range(int(steps)):
        elapsed = j*h+arb(0).union(h)
        residual = [v+arb(0, '.00003') for v in correction]
        wide_logs = [v+arb(0, '.001') for v in logs]
        inverse = [start+r+2*b*elapsed/LOOP for start, r, b in zip(inverse_high, residual, BETA)]
        assert min(inverse) > q('.8')
        gs = [1/x for x in inverse]
        ys = [x.exp() for x in wide_logs]
        tr = [sum(ys[i:i+3]) for i in (0, 3, 6)]
        total = 3*(tr[0]+tr[1])+tr[2]
        # Backward elapsed time, with the one-loop gauge term integrated exactly.
        dx = [2*sum(b*g-d*t for b, g, d, t in zip(br, gs, dr, tr))/LOOP**2
              for br, dr in zip(B, D)]
        dy = []
        for f, row in enumerate(C):
            other = ys[3*(1-f):3*(1-f)+3] if f < 2 else [arb(0)]*3
            gauge = sum(c*g for c, g in zip(row, gs))
            dy.extend(-(3*(u-v)+2*(total-gauge))/LOOP for u, v in zip(ys[3*f:3*f+3], other))
        assert all(h*abs(d) < q('.00003') for d in dx)
        assert all(h*abs(d) < q('.001') for d in dy)
        correction = [v+h*d for v, d in zip(correction, dx)]
        logs = [v+h*d for v, d in zip(logs, dy)]
    gauge = [(v+r+2*b*q(LOG_SCALE)/LOOP)**q('-1/2') for v, r, b in zip(inverse_high, correction, BETA)]
    yukawa = [(v/2).exp() for v in logs]
    gauge_centers = ('.467453', '.63811', '1.0583')
    yukawa_centers = ('.00000615', '.00311', '.8616', '.0000135', '.000268', '.01401', '.0000028227', '.00059465', '.010101')
    for values, centers, radius in [(gauge, gauge_centers, '.001'), (yukawa, yukawa_centers, '.01')]:
        print('Relative displacements of the 1 TeV inputs:', [v.str(4) for v in (x/q(c)-1 for x, c in zip(values, centers))])
        assert all(abs(v/q(c)-1) < q(radius) for v, c in zip(values, centers))


def shell_selection(s_upper):
    levels = {}
    for v in product(range(30), repeat=3):
        if sum(i*i for i in v) == 867:
            radius = sum((i-17)**2 for i in v)
            levels.setdefault(radius, []).append(v)
    counts = {radius: (len(pts), sum(x == y for x, y, z in pts)) for radius, pts in levels.items()}
    assert counts == {0: (1, 1), 68: (3, 1), 136: (9, 1), 544: (6, 0)}
    choices = {}
    for size in range(1, 5):
        for layers in combinations(sorted(counts), size):
            rank, trace = map(sum, zip(*(counts[r] for r in layers)))
            choices[layers] = F(trace, rank)
    assert len(choices) == 15 and choices[(136, 544)] == AMPLITUDE
    # Even projection of a second Weyl generation gives the same signed coefficient.
    assert F(8-7, 8+7) == AMPLITUDE

    def projection(step, inverse, derivative):
        if step >= 23*64:
            r1, r2 = inverse[0]/inverse[2], inverse[1]/inverse[2]
            reading = -r1+(5*r2-3)/2
            assert q('.18') < reading < q('.22') and 4*arb.pi()*inverse[2] > 35

    running['enclosure'](projection)
    # All signed choices, including unit distortion, share this coarse small-strength bound.
    x = 13*60*s_upper
    assert 4*x < 1
    mass = x/(1-x)
    moment = x/(1-4*x)
    vacuum = (1-mass)**3-mass**3
    assert vacuum > 0
    bare_lower = 4*s_upper/(1+60*s_upper)
    axis_error = moment*mass**2/(vacuum*bare_lower)
    ratio_error = (5+F('3.1'))*axis_error/(1-axis_error)
    assert 0 < axis_error < 1 and ratio_error < F('.00002')
    for layers, trace in choices.items():
        for t in (trace, -trace):
            if layers == (136, 544) and t == AMPLITUDE:
                continue
            separated = max(F('.18')-F('3.1')*t, F('3.1')*t-F('.22'))-ratio_error
            assert 35*separated/(5+F('3.1')*abs(t)+ratio_error) > F('.5')


def verify():
    backward_input()
    axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    rows, den = local['tables']((7, 13), (1+AMPLITUDE)/200)
    errors, _ = local['exterior'](rows, den, (7, 13, 13), axes, AMPLITUDE)
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    assert 4*(1+AMPLITUDE+max(errors)) < 5
    limits = {1: '1.094', 7: '.01561', 46: '.000133', 320: '.000000415',
              1000: '.0000000108', 3200: '.000000000226'}
    for index, limit in limits.items():
        cell = prefix[index-1][0]
        low = F(cell, 1 << 400)
        high = F(cell+1, 1 << 400)
        s_upper = 1/(F(1, 4)+low**2)
        assert 0 < s_upper < F(1, 200)
        c_target = 4*arb.pi()*q(COLOR)
        # C = 4 pi G3 / P with P=4s/[1+(30+2t)s].
        desired = c_target*(q(F(1, 4)+low**2).union(q(F(1, 4)+high**2))+q(30+2*AMPLITUDE))/4
        unit = 1 << 40
        floor = (desired*unit).floor().unique_fmpz()
        assert floor is not None
        calibration = F(int(floor), unit)
        rounding = q(F(5, 1)*s_upper/unit)
        error = c_target*q(max(errors)*(200*s_upper)**2)+rounding
        assert error < q(limit)
        if index == 3200:
            assert int(floor) == 158545267438980033232
            shell_selection(s_upper)
        print('Zero, dyadic calibration numerator over 2^40, maximum inverse-alpha error:',
              index, int(calibration*unit), error.str(9))


if __name__ == '__main__':
    verify()
