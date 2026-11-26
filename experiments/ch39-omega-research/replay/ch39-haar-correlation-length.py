"""Finite-volume conditional edge influences for the SU(5) weak-plane plaquette model."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from json import loads
from math import comb, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ZERO = (0, 0, 0, 0, 0)
T = F(1, 15)
KAPPA = F(1, 8)
STRIPS = [rows for k in range(1, 5) for rows in combinations(range(5), k)]


@cache
def tensor(weight):
    result = []
    for strip in STRIPS:
        p = [x+int(i in strip) for i, x in enumerate(weight)]
        if all(p[i] >= p[i+1] for i in range(4)):
            result.append(tuple(x-p[-1] for x in p))
    return result


def singlet_moments():
    states = {ZERO: 1}
    moments = []
    for degree in range(16):
        moments.append(states.get(ZERO, 0))
        if degree < 15:
            out = defaultdict(int)
            for weight, value in states.items():
                for child in tensor(weight):
                    out[child] += value
            states = out
    assert moments[:3] == [1, 0, 4]
    den = prod(j-i for i, j in combinations(range(5), 2))
    assert sum(m*prod(w[i]-w[j]+j-i for i, j in combinations(range(5), 2))
               for w, m in states.items()) == 30**15*den
    return moments


def influence(s, moments):
    h = 13
    num = [comb(h, k)*(s*(1+T))**k for k in range(h+1)]
    vacuum = 2-sum(a*b for a, b in zip(num, moments))
    assert vacuum > 0
    second = (1+KAPPA)/(1-KAPPA)*(1+T)**2*sum(a*b for a, b in zip(num, moments[2:]))/vacuum
    # Integer ceiling of the square root on the 2^-40 mesh.
    grid = 1 << 40
    rms = F(isqrt(second.numerator*grid**2//second.denominator)+1, grid)
    assert rms*rms > second and 1 > s*rms
    return 28*min((30+2*T)*s, s*rms/(1-s*rms))


def verify():
    moments = singlet_moments()
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    limits = {1: '.3893', 2: '.1623', 4: '.0751', 46: '.003735', 3200: '.000004866'}
    for index, limit in limits.items():
        low = F(prefix[index-1][0], 1 << 400)
        rate = influence(1/(low**2+F(1, 4)), moments)
        assert 0 < rate < F(limit) < 1
        if index == 3200:
            assert 3*rate**3 < F('3.46e-16')
        print('First allowed zero, conditional variation contraction:', index, float(rate))


if __name__ == '__main__':
    verify()
