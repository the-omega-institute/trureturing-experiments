"""External calibration versus integer microscopic amplification on the exact SU(5) center."""
from fractions import Fraction as F
from itertools import product
from json import loads
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
model = run_path(str(ROOT/'ch39-uniform-weak-center.py'))
q, arb = model['q'], model['arb']
T = model['AMPLITUDE']
C = F(158545267438980033232, 1 << 40)
POWER = round(C)


def verify():
    assert POWER == 144196081
    floor = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix'][3199][0]
    gamma = q(F(floor, 1 << 400)).union(q(F(floor+1, 1 << 400)))
    s = 1/(q(F(1, 4))+gamma**2)
    # cos(2 pi/5), cos(4 pi/5) in the real quadratic field.
    first, second = (arb(5).sqrt()-1)/4, (-arb(5).sqrt()-1)/4
    native = q(30+2*T)
    values = [2*(q(5-T)*a+q(10+2*T)*b) for a, b in [(first, second), (second, first)]]
    ratios = [(1+s*f)/(1+s*native) for f in values]
    r = F(1, 1 << 476)
    assert all(q('.99') < u < 1 and u**POWER < q(r) for u in ratios)
    assert ratios[1]**POWER > q(r/2)
    bare = 4*s/(1+s*native)
    high_color = 4*arb.pi()*q(model['COLOR'])
    mismatch = abs(POWER*bare-high_color)*q(1+T)
    assert mismatch < q('.000000084')
    # An anchored connected support of k faces has at most (4*39)^(k-1) embeddings.
    defect = 4*r/(1-624*r)
    assert defect < F(1, 1 << 473)
    assert 2*2**432*defect < F(1, 1 << 40)
    print('Microscopic copy count, bare inverse-alpha mismatch:', POWER, mismatch.str(9))
    print('Excited center face probability < 2^-473; for area <= 2^432, |1-<W>| < 2^-40.')

    # On a tetrahedron, enumerate the image of the coboundary and remove complete components.
    fluxes = set()
    for a, b, c in product(range(5), repeat=3):
        fluxes.add((a, b, c, (a-b+c) % 5))
    assert len(fluxes) == 125 and all(not any(f) or sum(x != 0 for x in f) >= 2 for f in fluxes)
    weights = [arb(1), ratios[0]**POWER, ratios[1]**POWER, ratios[1]**POWER, ratios[0]**POWER]
    from math import prod
    partition = sum(prod(weights[k] for k in f) for f in fluxes)
    excited = sum(prod(weights[k] for k in f) for f in fluxes if f[0])
    assert excited/partition < q(defect)


if __name__ == '__main__':
    verify()
