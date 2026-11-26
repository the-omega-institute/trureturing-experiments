"""A diagonal spectator U(1), its high-scale kinetic cost, and a relative-curvature completion."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import arb

ROOT = Path(__file__).resolve().parent
model = run_path(str(ROOT/'ch39-gauge-input-cover.py'))
running, q = model['running'], model['q']


def crossings(gauge):
    left = right = False
    allowed = []

    def observe(j, inv, deriv):
        nonlocal left, right
        if j < 23*64:
            return
        weak = inv[1]-q('16/15')*inv[2]
        left |= j == 26*64 and weak > 0
        right |= j == 28*64 and weak < 0
        assert deriv[1]-q('16/15')*deriv[2] < 0
        if 0 in weak:
            value = q('5/3')*(inv[0]-q('24/25')*inv[2])
            allowed.append((j, value))

    running['integrate'](gauge, 33*64, observe)
    assert left and right and allowed
    return allowed


def verify():
    centers = model['GAUGE']
    central = crossings([q(c) for c in centers])
    assert all(q('.01817') < value < q('.03425') and 26*64+45 <= j <= 26*64+48 for j, value in central)
    radius = F(1, 1000)
    for sides in product((0, 1), repeat=3):
        gauge = [q(c*(1+(side-1)*radius)).union(q(c*(1+side*radius))) for side, c in zip(sides, centers)]
        assert all(value < q('.07835') for j, value in crossings(gauge))
    assert 4*arb.pi()*q('.07835') < 1
    # Positive relative-curvature kinetic term k*(F_X - sqrt(3/5)*F_Y)^2.
    # On the unbroken graph F_X = sqrt(3/5)*F_Y; the extra energy vanishes identically.
    print('Central diagonal kinetic coefficient: 0.01817 < 1/g_X^2 < 0.03425.')
    print('Across the full input box: 1/g_X^2 < 0.07835, hence alpha_X > 1 if 1/g_X^2 > 0.')


if __name__ == '__main__':
    verify()
