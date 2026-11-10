"""Characteristic polynomial of the A233128 class transfer and its split against Hardin's polynomial.

Run `python3 hardin-a233128/enumerate.py` first; this reads its `data.json`.
"""
import json
from fractions import Fraction
from pathlib import Path

data = json.loads(Path(__file__).with_name('data.json').read_text())
T = data['T']
C = data['C']
n = len(T)

def charpoly(m):
    # Faddeev-LeVerrier over Q; returns monic coefficients, highest degree first.
    size = len(m)
    a = [[Fraction(x) for x in row] for row in m]
    mk = [[Fraction(0)] * size for _ in range(size)]
    coeffs = [Fraction(1)]
    for k in range(1, size + 1):
        prev = coeffs[-1]
        for i in range(size):
            mk[i][i] += prev
        mk = [[sum(a[i][t] * mk[t][j] for t in range(size)) for j in range(size)]
              for i in range(size)]
        coeffs.append(-sum(mk[i][i] for i in range(size)) / k)
    assert all(c.denominator == 1 for c in coeffs)
    return [int(c) for c in coeffs]

def divide(num, den):
    num = list(num)
    q = []
    while len(num) >= len(den):
        f = num[0] // den[0]
        assert f * den[0] == num[0]
        q.append(f)
        num = [x - f * y for x, y in zip(num, den + [0] * (len(num) - len(den)))][1:]
    return q, num

def value(p, x):
    v = 0
    for c in p:
        v = v * x + c
    return v

chi = charpoly(T)
hardin = [1] + [-c for c in C]
quotient, remainder = divide(chi, hardin)
assert all(r == 0 for r in remainder)
assert len(quotient) == 8 and value(quotient, 0) != 0
integer_roots = [s * 2 ** k for k in range(0, 60) for s in (1, -1) if value(quotient, s * 2 ** k) == 0]
# Integer roots of a monic integer polynomial divide its constant term, here -2^50.
assert abs(quotient[-1]) == 2 ** 50 and integer_roots == []
sums = [sum(v) for v in data['D']]
assert data['krylov_rank'] == 14
print('charpoly degree', len(chi) - 1)
print('quotient', quotient)
print('integer roots of quotient', integer_roots)

import itertools
classes_by_width = []
for w in range(2, 8):
    proper = [p for p in itertools.product(range(3), repeat=w) if all(x != y for x, y in zip(p, p[1:]))]
    assert len(proper) == 3 * 2 ** (w - 1)
    def key(p):
        b = tuple(int(p[j + 2] != p[j]) for j in range(w - 2))
        return min(b, b[::-1])
    perms = list(itertools.permutations(range(3)))
    orbits = {min(min(tuple(s[x] for x in q) for s in perms) for q in (p, p[::-1])) for p in proper}
    assert len(orbits) == len({key(p) for p in proper})
    classes_by_width.append(len(orbits))
assert classes_by_width == [1, 2, 3, 6, 10, 20]
print('relabelling-and-reversal classes, widths 2..7:', classes_by_width)
