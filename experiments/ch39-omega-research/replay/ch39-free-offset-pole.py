"""Exact interval checks for arbitrary six positive integer time offsets."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
previous = run_path(str(root / 'ch39-eight-time-pole-control.py'))
B, mul, rounded, sqrt_cell, node = (previous[k] for k in
    ('B', 'mul', 'rounded', 'root_cell', 'node'))
base = previous['base']
squares = tuple(j*j for j in range(1, 7))
sites = (1, 2, 3, 4, 11, 12, 3200)
z, witness = previous['quadrature'](list(sites), Q(0))
witness_weights = dict(zip(sites, witness))
top = B(1) / z[2]
scale = B(Q('0.000835'))


def preparation_bounds():
    numerator = base['data']('ch39-natural-mass-cell.json')
    mass = B(Q(numerator, 2**288), Q(numerator+1, 2**288))
    first = B(Q('.03276'), Q('.03294'))
    others = B(Q('.01565'))
    length = B(Q('.08063'))
    # The 20-coordinate weight vector is multilinear in the six parameters.
    # Concavity gives the common lower bound on the sum of square roots.
    for mask in range(64):
        lengths = []
        for i in range(1, 21):
            point = node(i)
            v = witness_weights.get(i, B(0))
            if mask and i != 3:
                factor = mul(*(point*(B(1)+B(b)*z[2])/(z[2]*(B(1)+B(b)*point))
                               for j, b in enumerate(squares) if mask & (1 << j)))
                natural = point / mass
                v = rounded(natural + factor * (v-natural))
            assert v.lo >= 0 and (i > 4 or v.lo > 0)
            lengths.append(sqrt_cell(point*v) if v.hi else B(0))
        assert first.lo < lengths[0].lo <= lengths[0].hi < first.hi
        assert max(v.hi for v in lengths[1:]) < others.lo
        assert sum(lengths, B(0)).lo > length.hi
    assert node(21).hi < others.lo**2
    assert first.lo**2 > scale.hi
    assert witness[0].hi < (z[0]/mass).lo
    assert (witness[0]*z[0]*previous['b41']).hi < scale.lo
    cosine = previous['cosines']
    loss = B(1) - sqrt_cell((B(1)+cosine[41])/B(2))
    assert length.lo > 2*first.hi + others.hi*loss.hi
    assert length.lo > 2*others.hi + first.hi*loss.hi
    sine = [B(0)] + [sqrt_cell((B(1)-v)/B(2)) for v in cosine[1:]]
    for j in range(1, 42):
        assert first.lo*(sine[j]-sine[j-1]).lo > others.hi*sine[1].hi


def crossing_bounds():
    model = previous['control']['cross']['ReadingBounds']()
    p, modes, q2, target = base['physical_inputs']()
    poles = [B(1)/node(i) for i in range(1, 5)]
    multiplicities = base['MODES']
    anchors = base['data']('ch39-free-offset-root-bounds.json')
    assert len(anchors) == 6
    upper_surface = B(650)
    higgs = B(125**2)*upper_surface/q2
    assert ((B(1)-p)/(B(1)-higgs*z[0])).hi < 1 and higgs.hi < poles[0].lo
    assert top.hi < 650 < poles[3].lo and (scale*upper_surface).hi < 8

    @lru_cache(None)
    def responses(y):
        natural = model.response(y)
        trial = sum((rounded(w/(B(1)-B(y)*v)) for v, w in zip(z, witness)), B(0))
        return natural, trial-natural

    def response(y, factors):
        natural, change = responses(y)
        return rounded(natural + factors*change)

    def secular(y, surface, mode, factors):
        return scale*(B(y)-B(125**2)*surface/q2)-mode+mode*(B(1)-p)*response(y, factors)

    for k, rows in enumerate(anchors):
        assert len(rows) == 3
        exposure = mul(*((top/t)**25 for t in poles[:2]))
        for j, row in enumerate(rows):
            for raw, mode, n in zip(row, modes, multiplicities, strict=True):
                y = Q(raw, 10**(5 if j else 12))
                assert (poles[j-1].hi if j else higgs.hi) < y < poles[j].lo
                assert (q2/upper_surface).lo*y > 1000**2
                # For k <= sum(lambda_j) <= k+1 bound the product using
                # concentration at endpoints below and arithmetic mean above.
                lower = 1+(top.lo-y)/(y+squares[-1])
                upper = 1+Q(k+1, 6)*(top.hi-y)/(y+squares[0])
                factors = B(lower**k, upper**6)
                assert secular(y, top, mode, factors).hi < 0
                exposure = mul(exposure, (top/B(y))**n)
        ratio = rounded(exposure**11/target)
        print('slab', k, 'upper ratio', float(ratio.hi))
        assert ratio.hi < Q('.399')

    # Upper bounds for the same eighteen roots at Y=650. Above pole three,
    # this surface also lies below the fourth root in each mixed block.
    upper_rows = [[2060, 1980, 1332, 2036, 1590, 2016],
                  [339, 336, 299, 338, 315, 337],
                  [569, 566, 526, 568, 546, 568]]
    exposure = mul(*((upper_surface/t)**25 for t in poles[:3]))
    for j, row in enumerate(upper_rows):
        for raw, mode, n in zip(row, modes, multiplicities, strict=True):
            y = Q(raw, 10**6 if j == 0 else 1)
            assert (poles[j-1].hi if j else 0) < y < poles[j].lo
            factor = mul(*((top+B(b))/B(y+b) for b in squares))
            assert secular(y, upper_surface, mode, B(1, factor.hi)).lo > 0
            exposure = mul(exposure, (upper_surface/B(y))**n)
    assert (exposure**11/target).lo > 1
    at_upper = mul(*((top+B(b))/B(650+b) for b in squares))
    for mode in modes:
        assert secular(Q(650), upper_surface, mode, B(at_upper.lo, 1)).hi < 0
    assert 3*25 + 3*sum(multiplicities) == 144
    print('Crossing between the third pole and 650, with exactly 144 open doublets.')


if __name__ == '__main__':
    preparation_bounds()
    crossing_bounds()
