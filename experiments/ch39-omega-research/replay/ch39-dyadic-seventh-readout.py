"""以有理数核验六次有限读数能确定的第七项二进制格。"""
import json
from fractions import Fraction as Q
from math import prod
from pathlib import Path


def multipliers(endpoint, numerators):
    nodes = [Q(v, 10**9) for v in numerators]
    slopes = [j*j-1 for j in range(2, 7)]

    def denominator(y):
        return prod(1+s*y for s in slopes)

    def kernel_numerator(y):
        return (1-y)*denominator(y)

    def nodal_product(y):
        return (y-endpoint)*prod((y-r)**2 for r in nodes)

    pole = -Q(1, 48)
    coefficient = kernel_numerator(pole)/nodal_product(pole)
    assert coefficient < 0

    def interpolant(y):
        return (kernel_numerator(y)-coefficient*nodal_product(y))/(1+48*y)

    constant = interpolant(Q(1))/denominator(Q(1))
    residues = []
    for slope in slopes:
        t = -Q(1, slope)
        residues.append(interpolant(t)/((1-t)*prod(1+s*t for s in slopes if s != slope)))
    return [constant, interpolant(Q(0))-constant-sum(residues), *residues]


def reading_bound(coeffs, inputs, upper):
    total = coeffs[0]
    for c, (lo, hi) in zip(coeffs[1:], inputs, strict=True):
        total += c*(hi if (c>0) == upper else lo)
    return total


def verify(file):
    data = json.loads(file.read_text())
    assert data['depth'] == 288
    scale = 2**data['depth']
    shift = Q('1.360900010e-8')
    intervals = []
    for j, n in enumerate(data['levels'][1:7], 1):
        delta = shift if j == 6 else 0
        intervals.append((Q(n, scale)+delta, Q(n+1, scale)+delta))
    lower = multipliers(1, (15, 1212594, 4902081))
    upper = multipliers(0, (1212564, 4902077, 997741029))
    bits = 44
    input_codes = (17563992027978, 17480697584871, 17345995515359,
                   17165516473595, 16946175576819, 16695471801235)
    cells = []
    for (a, b), n in zip(intervals, input_codes, strict=True):
        assert int(a*2**bits) == int(b*2**bits) == n
        cells.append((Q(n, 2**bits), Q(n+1, 2**bits)))
    low = reading_bound(lower, cells, False)
    high = reading_bound(upper, cells, True)
    assert 32071996606 < low*2**35 < high*2**35 < 32071996607
    reference_cell = data["levels"][7] // 2**(data["depth"]-35)
    assert reference_cell == 32071992476
    assert (data["levels"][7]+1) // 2**(data["depth"]-35) == reference_cell
    print("44位输入格", *input_codes)
    print("第七项35位整数与参考整数", int(low*2**35), reference_cell)


if __name__ == '__main__':
    import sys
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent
    verify(root/'ch39-minimum-tick-finite-inputs.json')
