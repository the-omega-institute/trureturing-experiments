"""Exact Taylor moments of symmetric real three-step symbols on the A4 torus."""

from collections import defaultdict
from itertools import combinations_with_replacement, permutations, product

import sympy as sp
from sympy.polys.matrices import DomainMatrix


FIELD = sp.QQ.algebraic_field(sp.sqrt(5))
C = FIELD.from_sympy((sp.sqrt(5)-1)/4)
COS = (FIELD.one, C, -C-FIELD.one/2, -C-FIELD.one/2, C)
SIN = (FIELD.zero, FIELD.one, 2*C, -2*C, -FIELD.one)
DIRECTIONS = ((1, -1, 0, 0, 0), (1, 1, -2, 0, 0), (1, 1, -1, -1, 0))


def representative(exponents):
    v = sorted(exponents)
    v = tuple(x-v[2] for x in v)
    return max(v, tuple(-x for x in reversed(v)))


def orbits(radius):
    return sorted({representative(v) for v in combinations_with_replacement(range(-radius, radius+1), 5)
                   if v[2] == 0 and sum(map(abs, v)) <= radius})


def moment(orbs, s, degree, direction):
    # The omitted derivative sign, factorial and (for odd orders) sin(2pi/5)
    # are common to the entire row, and do not affect equality or rank.
    row = []
    for v in orbs:
        orbit = set(permutations(v))
        value = sum(sum(x*y for x, y in zip(p, direction))**degree for p in orbit)
        phase = (COS if degree % 2 == 0 else SIN)[s*sum(v) % 5]
        row.append(phase*value/FIELD.convert(len(orbit)))
    return row


def rank(rows):
    return DomainMatrix(rows, (len(rows), len(rows[0])), FIELD).rank()


def differences(orbs, degrees):
    return [[x-y for x, y in zip(moment(orbs, 1, k, d), moment(orbs, 0, k, d))]
            for k in degrees for d in DIRECTIONS]


def verify():
    two, three = orbits(2), orbits(3)
    assert (len(two), len(three)) == (5, 9)
    constraints = differences(two, (0, 2, 3))
    heavy = [x-y for x,y in zip(moment(two, 2, 0, DIRECTIONS[0]), moment(two, 0, 0, DIRECTIONS[0]))]
    assert rank(constraints) == rank(constraints+[heavy]) == 3

    flat = [moment(three, s, k, d) for s in (0, 1) for k in (0, 2, 3, 4) for d in DIRECTIONS]
    assert rank(flat) == 9

    t = 2*(1+C)
    f = -(t+6*C)/(6*t+28*C)
    terms = [((), t-1), ((1,), 2*t*f), ((1, 1, -1), -t*(1+f)), ((2, -1), -t*f),
             ((1, 1), -5*f), ((2,), f), ((1, 1, 1), 1+4*f)]
    coefficients = defaultdict(lambda: FIELD.zero)
    for word, value in terms:
        for indices in product(range(5), repeat=len(word)):
            v = [0]*5
            for i, power in zip(indices, word):
                v[i] += power
            coefficients[representative(v)] += value/(5**len(word))
    assert set(coefficients).issubset(three)
    budget = sum(abs(FIELD.to_sympy(x)) for _, x in terms)
    assert sp.simplify(7-budget).is_positive

    def evaluate(row):
        return sum((a*coefficients[v] for a, v in zip(row, three)), FIELD.zero)

    assert all(evaluate(row) == FIELD.zero for row in differences(three, (0, 2, 3, 4)))
    centers = [evaluate(moment(three, j, 0, DIRECTIONS[0])) for j in range(3)]
    assert centers == [FIELD.zero, FIELD.zero, FIELD.from_sympy(5*(1+sp.sqrt(5))/4)]
    print("Range-two constraints and heavy-contrast augmented ranks: 3, 3.")
    print("Range-three absolute fourth-order flatness rank: 9 of 9.")
    print("Range-three selector matches light jets through order four, heavy value", FIELD.to_sympy(centers[2]))


if __name__ == "__main__":
    verify()
