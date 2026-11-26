"""Endpoint covering in the three initial gauge coordinates, uniformly in the Yukawa sector."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
witness = run_path(str(ROOT/'ch39-uniform-weak-center.py'))
running = witness['running']
q = running['q']
GAUGE = tuple(map(F, ('.467453', '.63811', '1.0583')))


def exclude_near_center(radius, lower):
    # In the common calibration certified for zero 3200, rescale by 4 pi G3.
    # Any other positive calibration is a scalar multiple of this same vector.
    color_inverse = 4*arb.pi()*q(witness['COLOR'])
    assert q('.000000000226')/color_inverse < q('1e-11')
    error = arb(0, '1e-11')
    v = [q(x)+error for x in (F(24, 25), F(16, 15), 1)]
    covered = []

    def observe(step, inverse, derivative):
        if step < 23*64:
            return
        gaps = [4*arb.pi()*abs(v[i]*inverse[j]-v[j]*inverse[i])/(v[i]+v[j])
                for i, j in combinations(range(3), 2)]
        assert any(x > q(lower) for x in gaps)
        covered.append(step)

    gauge = [q(c*(1-radius)).union(q(c*(1+radius))) for c in GAUGE]
    running['integrate'](gauge, 33*64, observe)
    assert covered == list(range(23*64, 33*64))
    print('Gauge-input half width, excluded inverse-alpha distance:', radius, lower)


def verify():
    t = witness['AMPLITUDE']
    target = [q(witness['COLOR']*x) for x in (1-F(3, 5)*t, 1+t, 1)]
    steps = witness['LOG_SCALE']*64
    assert steps.denominator == 1
    radius = F('.0007')
    gauge = [q(c*(1-radius)).union(q(c*(1+radius))) for c in GAUGE]
    for i, center in enumerate(GAUGE):
        for side in (-1, 1):
            face = gauge.copy()
            face[i] = q(center*(1+side*radius))
            _, endpoint = running['integrate'](face, int(steps))
            signed = -side*(endpoint[i]-target[i])
            print('Input-box face, inward inverse-gauge margin:', i+1, side, signed.str(8))
            assert signed > q('.00019')
    for radius, lower in [(F(0), '.0287'), (F('.0003'), '.0108'),
                          (F('.0004'), '.0066'), (F('.0005'), '.00017')]:
        exclude_near_center(radius, lower)


if __name__ == '__main__':
    verify()
