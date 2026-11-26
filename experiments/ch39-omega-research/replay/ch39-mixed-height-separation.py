"""Ordinary color-channel separation of independently assigned actual zero heights."""
from fractions import Fraction as F
from json import loads
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
witness = run_path(str(ROOT/'ch39-uniform-weak-center.py'))
q = witness['q']
t = witness['AMPLITUDE']
offset = 30+2*t


def verify():
    local = witness['local']
    rows, den = local['tables']((7, 13), (1+t)/200)
    errors, _ = local['exterior'](rows, den, (7, 13, 13), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), t)
    assert 0 <= min(errors) and max(errors) < F('.027')
    cells = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    assert len(cells) == 3200
    grid = 1 << 400
    heights = [(F(k[0], grid), F(k[0]+1, grid)) for k in cells]
    assert all(left[1] < right[0] for left, right in zip(heights, heights[1:]))
    lower = [1/(F(1, 4)+hi**2) for _, hi in heights]
    upper = [1/(F(1, 4)+lo**2) for lo, _ in heights]

    def bare(s):
        return 4*s/(1+offset*s)

    def relative(s):
        return F('.027')*(200*s)**2

    bounds = []
    for j in range(3199):
        ratio = bare(upper[j+1])/bare(lower[j])
        error = relative(upper[j])
        bounds.append(35*(1-ratio-2*error)/(1+ratio))
    worst = min(range(len(bounds)), key=bounds.__getitem__)
    assert worst+1 == 3150 and min(bounds) > F('.0014541')
    print('Smallest certified color gap:', worst+1, worst+2, float(bounds[worst]))

    # The lower envelope P(s)(1-a*s^2) increases throughout the allowed range.
    cap = F(1, 200)
    a = F('.027')*200**2
    assert 1-3*a*cap**2-2*offset*a*cap**3 > 0
    calibration = F(158545267438980033232, 1 << 40)
    color_target = 4*arb.pi()*q(witness['COLOR'])
    s = lower[-2]
    residual = q(calibration*bare(s)*(1-relative(s)))-color_target
    assert residual > q('.0197255')
    print('Fixed-calibration gap for any earlier zero:', residual.str(10))

    # Control with a dyadic strength that does not use a zero-height input.
    control = F(1, 1 << 24)
    unit = 1 << 40
    numerator = (unit*color_target/q(bare(control))).floor().unique_fmpz()
    assert numerator is not None
    residual = color_target*q(max(errors)*(200*control)**2)+q(5*control/unit)
    assert residual < q('.00000000016')
    print('Dyadic-strength control: common numerator over 2^40, error:', numerator, residual.str(9))


if __name__ == '__main__':
    verify()
