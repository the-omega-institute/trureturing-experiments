"""Constitutive metricity and distinct pair sums for the first 3200 actual zero stiffnesses."""
from fractions import Fraction as F
from heapq import heapify, heapreplace, heappop
from itertools import combinations, product
from json import loads
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GRID = 1 << 112


def metric_checks():
    pairs = list(combinations(range(4), 2))
    shapes = {}
    wedge = sp.zeros(6)
    for i, a in enumerate(pairs):
        for j, b in enumerate(pairs):
            if not set(a).intersection(b):
                order = a+b
                wedge[i, j] = (-1)**sum(x > y for x, y in combinations(order, 2))
    for labels in product(range(3), repeat=4):
        if 1 not in labels or 2 not in labels:
            continue
        bivector = sp.Matrix([int((labels[i],labels[j]) == (1,2))-int((labels[i],labels[j]) == (2,1))
                              for i, j in pairs])
        key = tuple(sorted((labels.count(1), labels.count(2))))
        shapes[key] = shapes.get(key, sp.zeros(6))+bivector*bivector.T/4
    full = sum(shapes.values(), sp.zeros(6))
    parallel = (full-sp.Rational(3, 2)*sp.eye(6))/6
    perpendicular = sp.eye(6)-parallel
    assert parallel**2 == parallel and parallel.rank() == 3
    weights = sp.symbols('a b c d')
    combined = sum((w*shapes[k] for w,k in zip(weights, sorted(shapes))), sp.zeros(6))
    a, b, c, d = weights
    assert combined == (a/2+b)*perpendicular+(a/2+3*b+2*c+2*d)*parallel
    z, lam = sp.symbols('z lam')
    delta = sp.zeros(6)
    delta[0, 0] = (z-1)/2
    closure = (wedge*(full+delta))**2
    baseline = sp.Rational(45, 4)
    assert sp.factor(closure.charpoly(lam).as_expr()-(lam-baseline)**4*(lam-baseline*(4+z)/5)**2) == 0


def verify():
    metric_checks()
    prefix = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    offset = F(1, 4)+30+F(2, 15)
    weights = []
    for cell, *_ in reversed(prefix):
        bottom, top = (F(x, 1 << 400) for x in (cell, cell+1))
        low, high = 4/(offset+top**2), 4/(offset+bottom**2)
        code = low*GRID//1
        assert high*GRID//1 == code
        weights.append(code)
    n = len(weights)
    assert n == 3200 and all(a < b for a, b in zip(weights, weights[1:]))
    heap = [(2*v, i, i) for i, v in enumerate(weights)]
    heapify(heap)
    previous = None
    best = None
    closest = None
    second_best = None
    count = 0
    while heap:
        total, i, j = heap[0]
        if previous is not None:
            gap = total-previous[0]
            if best is None or gap < best:
                second_best = best
                best = gap
                closest = (previous[1], tuple(sorted((n-i, n-j))))
            elif second_best is None or gap < second_best:
                second_best = gap
        previous = (total, tuple(sorted((n-i, n-j))))
        count += 1
        if j+1 < n:
            heapreplace(heap, (weights[i]+weights[j+1], i, j+1))
        else:
            heappop(heap)
    assert count == n*(n+1)//2 and best > 2
    assert second_best > best+4 and closest == ((2112, 2576), (1977, 2863))
    separation = F(best-2, GRID)
    assert separation > F('3.53258e-19')
    # In either direction, a single anomalous coordinate square splits the closure eigenvalues.
    fraction = min(F(b, a+1)-1 for a, b in zip(weights, weights[1:]))
    split = fraction/(5+4*fraction)
    assert split > F('1.66186e-5')
    print('Pair sums, certified separation and closest grid pair:', count, float(separation), closest)
    print('Smallest certified relative closure splitting:', float(split))


if __name__ == '__main__':
    verify()
