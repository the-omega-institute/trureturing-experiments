"""Gauge-factor comparisons with common positive stiffness shares on symmetric powers."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
from flint import arb, ctx

root = Path(__file__).resolve().parent
base = run_path(str(root / 'ch39-native-arc-gauge-blocking.py'))
irrep = base['irrep']


def ball(q):
    return arb(q.numerator)/q.denominator


def examples():
    # Height enclosures are dyadic records from the same archived zero prefix.
    import json
    prefix = json.loads((root / 'ch39-tick-separators.json').read_text())['prefix']
    ctx.prec = 256
    raw = []
    for cell, *_ in prefix[:7]:
        height = arb(cell).union(arb(cell+1))/2**400
        raw.append(1/(arb(1)/4 + height**2))
    shares = [w/sum(raw) for w in raw]
    values = []
    for group, k, unit in [('U1', Q(1,4), Q(1,60)), ('SU2', Q(1,60), Q(1,4)),
                           ('SU3', Q(1,40), Q(1,6))]:
        terms = []
        for n, w in enumerate(shares, 1):
            c, dim = irrep(group, n)
            terms.append((ball(c), dim, ball(k)*w/ball(c)))
        vacuum = 1 - sum(p for _, _, p in terms)
        den = vacuum**2 + sum(p**2 / dim for c, dim, p in terms)
        values.append(ball(unit)*sum(c*p**2/dim for c, dim, p in terms)/den)
    print('First seven RH weights, two-cell normalized stiffness', [a.str(13) for a in values])
    assert values[0] > 30*values[1] and values[1] > ball(Q(164274,55700))*values[2]


def comparisons():
    # Uniform constants for 0 < kappa < 1/240 and all integers b >= 2.
    assert Q(3,2)*(Q(3,5)-Q(1,342)) > Q(30,29)*(Q(25,24)*Q(3,5)+Q(1,116))
    assert Q(3,2) > Q(30,29)*Q(25,24)
    for cells in (2, 3, 8):
        z2 = (1-Q(1,60))**cells + Q(1,60)**cells / 4**(cells-1)
        z3 = (1-Q(1,40))**cells + Q(1,40)**cells / 18**(cells-1)
        if cells == 2:
            assert 3*z3/z2 == Q(164274,55700)
        assert 3**(cells-1)*z3/z2 > 1


if __name__ == '__main__':
    comparisons()
    examples()
