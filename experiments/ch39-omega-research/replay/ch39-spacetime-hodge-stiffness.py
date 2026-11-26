"""Constant simplex fluxes, periodic correctors, and a common rational metric in 4D."""
from collections import defaultdict
from itertools import combinations, product
from math import factorial
import sympy as sp


def check(d, p):
    axes = list(combinations(range(d), p))
    gram = sp.zeros(len(axes))
    codifferential = defaultdict(lambda: sp.zeros(len(axes), 1))
    total = 0
    for labels in product(range(p+1), repeat=d):
        if not all(k in labels for k in range(1, p+1)):
            continue
        blocks = tuple(tuple(i for i, c in enumerate(labels) if c == k) for k in range(1, p+1))
        signs = []
        for face in axes:
            values = [labels[i] for i in face]
            sign = 0
            if sorted(values) == list(range(1, p+1)):
                sign = (-1)**sum(x > y for x, y in combinations(values, 2))
            signs.append(sign)
        flux = sp.Matrix(signs)
        gram += flux*flux.T
        codifferential[blocks[1:]] += flux
        codifferential[blocks[:-1]] += (-1)**p*flux
        for j in range(p-1):
            merged = tuple(sorted(blocks[j]+blocks[j+1]))
            boundary = blocks[:j]+(merged,)+blocks[j+2:]
            codifferential[boundary] += (-1)**(j+1)*flux
        total += 1
    assert all(v.is_zero_matrix for v in codifferential.values())
    tensor = sp.eye(d)+sp.ones(d)
    compound = sp.Matrix([[tensor.extract(a, b).det() for b in axes] for a in axes])
    gram /= factorial(p)**2
    assert gram == sp.Integer(p+1)**(d-p-1)*compound/factorial(p)
    return gram, total


def verify():
    for d in range(2, 6):
        for p in range(1, d+1):
            gram, total = check(d, p)
            if d == 4:
                print('Form degree, simplex orbits, eigenvalues:', p, total, gram.eigenvals())
            if (d, p) == (4, 2):
                assert total == 50 and gram.eigenvals() == {sp.Rational(3, 2): 3, sp.Rational(15, 2): 3}
    metric = sp.eye(4)-sp.ones(4)/5
    assert metric.inv() == sp.eye(4)+sp.ones(4) and metric.det() == sp.Rational(1, 5)
    steps = [sp.eye(4)[:, i] for i in range(4)]+[-sp.ones(4, 1)]
    five = sp.Matrix([[(x.T*metric*y)[0] for y in steps] for x in steps])
    assert five == sp.eye(5)-sp.ones(5)/5
    # Four-dimensional shape classes: off-diagonal terms cannot cancel with positive weights.
    shape = defaultdict(int)
    for lab in product(range(3), repeat=4):
        if 1 in lab and 2 in lab:
            def coord(a, b):
                return int((lab[a], lab[b]) == (1, 2))-int((lab[a], lab[b]) == (2, 1))
            shape[tuple(sorted((lab.count(1), lab.count(2))))] += coord(0, 1)*coord(0, 2)
    assert shape == {(1, 1): 0, (1, 2): 2, (1, 3): 2, (2, 2): 2}


if __name__ == '__main__':
    verify()
