"""Arithmetic of independently labelled Weyl faces and the one-loop running plane."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, lcm, prod
from pathlib import Path
from runpy import run_path
import json

ROOT = Path(__file__).resolve().parent
ORDER = 802
MODULUS = 1 << 40


def scale_normal():
    em, sw = F('127.930'), F('0.23129')
    start = (3*(1-sw)*em/5, sw*em, 1/F('0.1180'))
    beta = (F(41, 10), F(-19, 6), F(-7))
    axis = [start[(j+1)%3]*beta[(j+2)%3]-start[(j+2)%3]*beta[(j+1)%3] for j in range(3)]
    common = lcm(*(x.denominator for x in axis))
    ints = [int(x*common) for x in axis]
    normal = tuple(a//gcd(*ints) for a in ints)
    assert normal == (-1595534194915, 3962837583051, -2727239411259)
    assert sum(normal) == -359936023123 and all(a % 2 for a in normal)
    for v in (start, beta):
        assert sum(a*b for a, b in zip(normal, v, strict=True)) == 0
    return normal


def unique_bits(rows, bits, collision):
    assert len({v % (1 << bits) for _, v in rows}) == len(rows)
    values = dict(rows)
    assert values[collision[0]] % (1 << (bits-1)) == values[collision[1]] % (1 << (bits-1))


def verify():
    prefix = json.loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    denoms = [(1 << 800)+(2*cell+1)**2 for cell, *_ in prefix]
    units = [pow(d, -1, MODULUS) for d in denoms]
    assert len(units) == 3200 and all(u % 8 == 1 for u in units)
    unique_bits([(j+1, u*u % MODULUS) for j, u in enumerate(units)], 28, (1406, 2626))
    for n, bits, pair in [(7, 18, ((1, 6, 6), (4, 4, 5))),
                          (15, 23, ((2, 3, 3), (8, 9, 15))),
                          (46, 33, ((5, 23, 26), (15, 28, 29)))]:
        rows = [(ids, prod(units[i-1] for i in ids) % MODULUS)
                for ids in combinations_with_replacement(range(1, n+1), 3)]
        unique_bits(rows, bits, pair)
        print(n, 'eligible labels;', len(rows), 'triples;', bits, 'low bits of the normalized residue.')
    # Direct heterogeneous tetrahedron check with three different labels and a separate base label.
    witnesses = [F(1 << ORDER, d) for d in denoms[:4]]
    cap = run_path(str(ROOT/'ch39-weyl-surface-fillings.py'))['cap'](prod(witnesses[1:]))
    base = 4*witnesses[0]/(1+30*witnesses[0])
    scaled = 18*base*((base+cap[1])/(base+cap[2])-1)/2**(3*ORDER)
    diff = scaled-prod(F(1, d) for d in denoms[1:4])
    val2 = run_path(str(ROOT/'ch39-weyl-uniform-crossing.py'))['valuation2']
    assert val2(diff) >= 2*ORDER-46
    ratio23 = (base+cap[1])/(base+cap[2])
    assert val2(ratio23-1) == 1601
    scaled_ratio = (ratio23-1)/2**1601
    assert (scaled_ratio.numerator*pow(scaled_ratio.denominator, -1, 8)) % 8 == 1
    normal = scale_normal()
    ratios = [(base+k)/(base+cap[2]) for k in cap]
    assert val2(sum(a*b for a, b in zip(normal, ratios))) == 0
    print('All one-loop scales excluded. Sharp 28 bits for 3200 homogeneous labels; cap-sum precision', 2*ORDER-46)


if __name__ == '__main__':
    verify()
