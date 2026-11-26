"""A phase strip with a fixed Wilson coefficient and a dynamically selected parallel axis."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
unit = runpy.run_path(str(root / 'ch39-unit-pentagon.py'))
m = unit['model']
q = m['real']


def weighted_terms():
    degree = 20
    first = Fraction(1, 2**10)
    cutoff = 2**13
    nodes = [arb.legendre_p_root(degree, j, weight=True) for j in range(degree)]
    yield q(first**2/2), m['schur_terms'](0, first)
    for power in range(-10, 13):
        half = q(Fraction(2)**power) / 2
        for rad, weight in nodes:
            k = half * (3 + rad)
            yield half * weight * k, m['schur_terms'](*unit['node_enclosure'](k))
    # Uniform in c including both collinear endpoints. The beta second derivative
    # bound gives |k F'_k| < 6000 on every complex Taylor disk used here.
    error = q(4 * 6000 * cutoff) / 2**(2 * degree) + m['tail_bound'](cutoff)
    yield arb(0, error.upper()), None


def gradient(terms, u, v, c, wilson_weight=1):
    su, x = u.sin_cos()
    sv, y = v.sin_cos()
    s = su * sv
    z = x*y + c*s
    # Expand the two-dimensional Schur determinant, differentiating only z.
    factors = (1-x*x, 1-y*y, 1-z*z, 2*(z-x*y), 2*(x*z-y), 2*(x-y*z))
    total = q(0)
    for weight, sectors in terms:
        if sectors is None:
            total += weight
            continue
        summed = q(0)
        for const, hs in sectors:
            a, b, d = (h.imag for h in hs)
            tr = const + 2*(hs[0].real*x + hs[1].real*y + hs[2].real*z)
            norm = sum(w*t for w,t in zip((a*a,b*b,d*d,a*b,a*d,b*d),factors))
            det = tr*tr - 4*norm
            assert det > 0
            summed += (4*tr*hs[2].real - 8*(-d*d*z+a*b+a*d*x-b*d*y))/det
        total += weight * s * summed
    scale = unit['pentagon']()[1] / (m['su']*m['sv'])**2
    return total + wilson_weight * scale * s * z


def strip_monotonicity():
    terms = list(weighted_terms())
    u = q('2.8') + arb(0, q('.0001'))
    worst = q(-10)
    for j in range(32):
        v = q('.3') + q(Fraction(j,64)).union(q(Fraction(j+1,64)))
        for i in range(8):
            c = q(Fraction(i-4,4)).union(q(Fraction(i-3,4)))
            val = gradient(terms, u, v, c)
            worst = max(worst, val.upper())
    print('strip upper total slope', worst.str(12), flush=True)
    assert worst < q('-.3389')


def gauge_faces():
    f = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
    tau0, tau1 = q('26.97'), q('27.03')
    tube = tau0.union(tau1)
    u = q('2.8')
    for v, sign in ((q('.3'), 1), (q('.8'), -1)):
        xs = f['spectrum'](u, q('.0001'), True, v, q(0))
        local = tube
        for j in range(3):
            t, res, _, _ = f['run'](local, xs, 8192 if j<2 else 16384, return_target=True)
            assert t.overlaps(local)
            local = local.intersection(t)
        print('gauge phase face margin', (sign*res).lower().str(12), flush=True)
        assert sign*res > 0
    lo = q(1)
    hi = q(1)
    for j in range(16):
        halfwidth = q(1)/64
        v = q('.3') + (2*j+1)*halfwidth
        xs = f['spectrum'](u, halfwidth + q('.0001'), False, v, q(0))
        lo = min(lo, f['run'](tau0, xs, 1024)[0].lower())
        hi = min(hi, (-f['run'](tau1, xs, 1024)[0]).lower())
    print('time face margins', lo.str(12), hi.str(12), flush=True)
    assert lo>0 and hi>0


if __name__ == '__main__':
    strip_monotonicity()
    gauge_faces()
