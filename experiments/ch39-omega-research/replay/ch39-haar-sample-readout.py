"""Population conjugation contrasts and finite-sample bounds for hidden face strength."""
from fractions import Fraction as F
from itertools import combinations
from json import loads
from pathlib import Path

ROOT = Path(__file__).resolve().parent
T = F(1, 15)
SUP = 30+2*T
GRID = 1 << 100


def down(x):
    return F((x*GRID)//1, GRID)


def up(x):
    return -down(-x)


def verify():
    z = (1, 1, 1, -1, -1)
    other = (-1, -1, 1, 1, 1)
    wedge = lambda a: tuple(a[i]*a[j] for i, j in combinations(range(5), 2))
    dims = (5, 10)
    tables = tuple(zip((z, wedge(z)), (other, wedge(other)), dims))
    haar_delta = 2*T*T*sum(F(sum((x-y)**2 for x, y in zip(a, b)), d) for a, b, d in tables)
    haar_native = 2*sum(sum((1-T*x)**2 for x in a)/d for a, _, d in tables)
    assert haar_delta == F(16, 375) and haar_native == 4*(1+T*T)
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    strengths = [(down(1/(F(1, 4)+F(cell+1, 1 << 400)**2)),
                  up(1/(F(1, 4)+F(cell, 1 << 400)**2))) for cell, *_ in prefix]
    mids = [(left[0]+right[1])/2 for left, right in zip(strengths, strengths[1:])]
    assert all(low > high for (low, _), (_, high) in zip(strengths, strengths[1:]))
    sufficient = []
    for i, (lo, hi) in enumerate(strengths):
        # Extend the endpoint labels with an artificial outer test to force a decreasing empirical line.
        larger = mids[i-1] if i else hi+(lo-mids[0])
        smaller = mids[i] if i < len(mids) else lo-(mids[-1]-hi)
        assert 0 < smaller < lo <= hi < larger < F(1, 100)
        margin = min(larger-hi, lo-smaller)
        b = SUP*hi
        weight_ratio = up(((1+b)/(1-b))**33)
        j_lower = haar_delta/(weight_ratio*(1+b))
        mean_gap = margin*j_lower/2
        # E_p[A_a^2] for A_a=(f-f')*(1+a*f'), followed by Chebyshev on two tests.
        variance = weight_ratio*(1+b)/(1-b)*(1+larger*SUP)**2*haar_delta
        sufficient.append(40*variance/mean_gap**2)
    hard = max(range(len(sufficient)), key=sufficient.__getitem__)

    a, b = 3149, 3150
    delta = strengths[a][1]-strengths[b][0]
    chi2 = delta**2*haar_native/(1-SUP*strengths[b][1])
    necessary = 2*F('.9')**2/chi2
    assert hard+1 == 3150 and max(sufficient) < 4*10**26
    assert necessary > F('1.0738e22')
    print('Sufficient independent draws, worst label:', float(max(sufficient)), hard+1)
    print('Necessary independent draws for 95% on both 3150 and 3151:', float(necessary))


if __name__ == '__main__':
    verify()
