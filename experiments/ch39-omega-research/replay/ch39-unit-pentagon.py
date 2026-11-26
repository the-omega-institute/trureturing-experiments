"""Unit Wilson coupling on the continuous shell, using Gauss integration with a Cauchy bound."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import acb, arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / 'ch39-gaussian-momentum.py'))
q = model['real']
shell = model['shell']


def pentagon():
    phases = {(i, j): z for i, j, _, z in shell['edges']}
    product = acb(1)
    for i, j in ((0, 7), (7, 4), (4, 6), (6, 8), (8, 0)):
        product *= phases[i, j] if i < j else phases[j, i].conjugate()
    s = model['su'] * model['sv']
    return -model['x'] * model['y'] / s, -4 * product.real * s**2


def node_enclosure(node):
    # Each Gauss node is enclosed between exact dyadic rationals for the edge coefficient cache.
    scale = 2**100
    low = Fraction(int((node * scale).lower().floor().unique_fmpz()), scale)
    high = Fraction(int((node * scale).upper().ceil().unique_fmpz()), scale)
    return low, high


def slopes_at(cosines):
    degree = 26
    epsilon = Fraction(1, 2**14)
    cutoff = 2**20
    low_terms = model['schur_terms'](0, epsilon)
    totals = [model['log_derivative'](low_terms, c) * q(epsilon**2 / 2) for c in cosines]
    assert all(abs(c) < q('.5') for c in cosines)
    nodes = [arb.legendre_p_root(degree, j, weight=True) for j in range(degree)]
    for power in range(-14, 20):
        left = q(Fraction(2)**power)
        halfwidth = left / 2
        midpoint = 3 * halfwidth
        for rad, weight in nodes:
            node = midpoint + rad * halfwidth
            terms = model['schur_terms'](*node_enclosure(node))
            for j, c in enumerate(cosines):
                totals[j] += halfwidth * weight * node * model['log_derivative'](terms, c)
    # On |k-3a/2|<=a, accretivity and the six moving edges give |k F'_k|<150.
    # Positive Gauss weights and exactness through 2n-1 bound both truncation errors.
    error = 4 * q(150) * cutoff / 2**(2 * degree) + model['tail_bound'](cutoff)
    return [total + arb(0, error.upper()) for total in totals]


if __name__ == '__main__':
    center, curvature = pentagon()
    assert curvature - q('.19632') > q('2.495')
    cosines = [q(x) for x in ('-.1567', '-.1559')]
    slopes = slopes_at(cosines)
    derivatives = [slope + curvature * (c - center) for c, slope in zip(cosines, slopes)]
    assert derivatives[0] < q('-.00013') and derivatives[1] > q('.000036')
    print('total slopes at faces', derivatives[0].upper().str(12), derivatives[1].lower().str(12))
    gauge = runpy.run_path(str(root / 'ch39-eight-zero-short-arc.py'))
    residuals = []
    times = []
    for j in range(8):
        left = (cosines[0] + q(j)/10000).acos()
        right = (cosines[0] + q(j+1)/10000).acos()
        angle = (left + right)/2
        width = (left - right)/2 + q('.00001')
        tau, residual = gauge['project_fixed_axis'](angle, width)
        assert residual < 0
        times.append(tau)
        residuals.append(residual)
    print('all candidate times', min(t.lower() for t in times).str(12), max(t.upper() for t in times).str(12))
    print('all candidate gauge residuals', min(r.lower() for r in residuals).str(12),
          max(r.upper() for r in residuals).str(12))
