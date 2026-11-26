"""A weak-center deformation matching two gauge ratios for a fixed measured triangle."""
from fractions import Fraction as F
from json import loads
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
local = run_path(str(ROOT/'ch39-exponential-face-stiffness.py'))
running = run_path(str(ROOT/'ch39-two-loop-shell-separation.py'))
q = running['q']
DELTA = F('.085')
AXES = ((5, 3, -8), (-2, 5, -3), (0, 0, 1))


def running_rectangle():
    samples = []

    def observe(step, inverse, derivative):
        if not 26*64 <= step <= 28*64:
            return
        r1, r2 = inverse[0]/inverse[2], inverse[1]/inverse[2]
        for i, ratio in enumerate((r1, r2)):
            assert abs((derivative[i]-ratio*derivative[2])/inverse[2]) < 1
        a = 5*r1+3*r2-8
        b = -r1+(5*r2-3)/2
        samples.append(step)
        if step == 26*64:
            assert a > q('.189')
        if step == 28*64:
            assert a < q('-.27')
        assert q('.191') < b < q('.213')

    running['enclosure'](observe)
    assert samples == list(range(26*64, 28*64+1))


def verify():
    running_rectangle()
    signs = (1, -1, -1, 1, -1)
    indices = local['INDICES']
    assert tuple(map(sum, zip(*indices))) == (20, 20, 20)
    assert tuple(sum(weight*v for weight, v in zip(signs, col)) for col in zip(*indices)) == (-12, 20, 0)
    variation = (F(-3, 5), F(1), F(0))
    assert sum(a*v for a, v in zip(AXES[0], variation)) == 0
    slope = sum(a*v for a, v in zip(AXES[1], variation))/2
    assert slope == F('3.1')
    rows, den = local['tables']((7, 13), (1+DELTA)/200)
    errors, _ = local['exterior'](rows, den, (7, 13, 13), AXES, DELTA)
    # The actual positive zero is enclosed in its 400-bit cell.
    ea, e2b, ec = errors
    # Uniform throughout [0, DELTA], including the endpoints needed for Markov's inequality.
    assert ec < F('.007')
    assert ea/(1-ec) < F('.154')
    assert slope*DELTA+(e2b/2+slope*DELTA*ec)/(1-ec) < F('.312')
    # Inverting the two projective coordinates then bounds both color ratios by two.
    assert (5*F('.154')+6*F('.312'))/31 < 1
    assert (2*F('.154')+10*F('.312'))/31 < 1
    assert 1+30*(1+DELTA)/200 < F(6, 5)
    assert 1-30*(1+DELTA)/200 > F(4, 5)
    assert 4*(1-ec)/F(6, 5) > 3 and 4*(1+ec) < 5
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    for index, lower, upper in [(1, '.046', '.085'), (2, '.058', '.072'), (7, '.061', '.069')]:
        gamma = F(prefix[index-1][0], 1 << 400)
        strength_cap = 1/(F(1, 4)+gamma*gamma)
        assert strength_cap < F(1, 200)
        factor = (200*strength_cap)**2
        ea, e2b, ec = [factor*e for e in errors]
        assert 0 <= ec < 1
        a_error = ea/(1-ec)
        b_error = (e2b/2+slope*DELTA*ec)/(1-ec)
        assert a_error < F('.189')
        assert slope*F(lower)+b_error < F('.191')
        assert slope*F(upper)-b_error > F('.213')
        print('First allowed zero, certified open parameter interval, A/B errors:',
              index, (lower, upper), float(a_error), float(b_error))


if __name__ == '__main__':
    verify()
