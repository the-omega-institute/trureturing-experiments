"""Positive uniform face weights and geometric multiplicities near the book's crossing."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path
import json

ROOT = Path(__file__).resolve().parent
patch = run_path(str(ROOT/'ch39-weyl-surface-fillings.py'))


def orient(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def valuation2(x):
    def order(n):
        n = abs(n)
        return (n & -n).bit_length()-1
    assert x
    return order(x.numerator)-order(x.denominator)


def filled_congruence(t):
    k = valuation2(t)
    assert k == 802 and 4*(k-9)-2 > 3*k-1
    target = patch['TARGET']
    assert valuation2(target-1) == 0
    bare = 4*t/(1+30*t)
    result = [bare+value for value in patch['cap'](t**3)]
    for i, j in ((0, 2), (1, 2)):
        assert valuation2(result[i]/result[j]-1) == 2*k-3 == 1601
    assert valuation2(result[0]/result[1]-1) == 2*k+1 == 1605
    print('Filled triangle: sharp universal exponent at the archived input is', 2*k-3)


def verify():
    cell = json.loads((ROOT/'ch39-tick-separators.json').read_text())['prefix'][0][0]
    h = F(2*cell+1, 2**401)
    t = 1/(F(1, 4)+h*h)
    assert F(1, 201) < t < F(1, 200)
    filled_congruence(t)
    _, weak, color = patch['crosscap'](t, t)
    _, t0, _ = patch['match'](t)
    common = 4*t0/(1+30*t0)
    bare = 4*t/(1+30*t)
    bands = round(bare/common)
    z = bands*color/(F(13, 10)-F(14, 3)*bands*color)
    third_face = z/t**2
    boundary_face = bands*common/(4-30*bands*common)
    cap = patch['cap'](z)
    result = [x+bands*m+4*boundary_face/(1+30*boundary_face)
              for x, m in zip(cap, (0, weak, color), strict=True)]
    assert result[0] == result[2] and result[1] == patch['TARGET']*result[0]
    assert F('.82885') < third_face/t < F('.82886')
    assert abs(boundary_face/t-1) < F(1, 10**86)
    print('Bands with two moderate fitted strengths:', bands)

    # Box about an identity translate in the image of the 15-dimensional ledger.
    delta = F(1, 10**9)
    assert F(450, 2*2**40) < delta/2
    assert 1+F(14, 3)*t**3 < 2 and (1+30*t)/(1-30*t) < 16
    block = 2*10**9*201**3
    assert block*t**3*delta/2 > 1
    reps = 9500*block
    assert 16*F(2)**(21*450-9500) < delta
    assert 6*delta+2*delta < F(1, 10**8)
    tiny = F(1, 10**8)
    box = F(1, 10**7)
    assert 1-box < (1-tiny)**2 and (1-tiny)**(-2) < 1+box
    disk = (F(142, 25), F(6, 5))
    band = (F(0), weak/color)
    assert band[1] == F(134217728000000008, 11438396227485)
    for near in product((1-box, 1+box), repeat=2):
        for target in product((1-F(1, 10**6), 1+F(1, 10**6)), (F('1.08'), F('1.10'))):
            assert orient(disk, band, target) > 0
            assert orient(band, near, target) > 0
            assert orient(near, disk, target) > 0
    print('Cap copies per concentrated face (sufficient):', reps)
    print('All rational targets in [1-1e-6,1+1e-6] x [1.08,1.10] lie strictly inside the three-patch triangle.')


if __name__ == '__main__':
    verify()
