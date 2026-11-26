"""Exact checks for the nine-point Behrend shell and a commuting rational clock.

For the geometry, run this file with Python and SymPy. For the arithmetic
origin of the dyadic clock table run ch39-five-tick-record-origin.py
(python-flint). ShellProjection.lean replays the finite integer compressions
and one oriented triple with Lean 4.35.0 and only Init.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


def multiply_by_vector(vector, identity, *, right):
    def column(q):
        tail = q[1:4, 0]
        cross = tail.cross(vector) if right else vector.cross(tail)
        return sp.Matrix.vstack(sp.Matrix([-vector.dot(tail)]), q[0] * vector + cross)

    return sp.Matrix.hstack(*(column(identity[:, j]) for j in range(4)))


def norm_obstruction_at(value, prime):
    """A prime congruent to 3 mod 4 occurs evenly in every nonzero rational sum of two squares."""
    assert prime % 4 == 3 and sp.isprime(prime)
    valuation = sp.factorint(sp.numer(value)).get(prime, 0) - sp.factorint(sp.denom(value)).get(prime, 0)
    return valuation % 2 == 1


def main():
    identity = sp.eye(4)
    points = [sp.Matrix(v) for v in product(range(5), repeat=3) if sum(x*x for x in v) == 17]
    assert len(points) == 9
    right = multiply_by_vector(points[0], identity, right=True)
    left = [multiply_by_vector(p, identity, right=False) for p in points]
    planes = [(identity + op*right/17)/2 for op in left]
    assert right**2 == -17*identity
    assert all(p == p.T and p*p == p and p.trace() == 2 for p in planes)
    for i, j in product(range(9), repeat=2):
        k = (17 + points[i].dot(points[j]))/34
        assert planes[i]*planes[j]*planes[i] == k*planes[i]
    for i, j, k in permutations(range(9), 3):
        p, q, r = (points[t] for t in (i, j, k))
        even = (17 + p.dot(q) + q.dot(r) + r.dot(p))/68
        odd = sp.Matrix.hstack(p, q, r).det()/1156
        loop = planes[i]*planes[j]*planes[k]*planes[i]
        assert loop == (even*identity + odd*right)*planes[i]

    shortest = sp.Rational(33, 34)**3
    complement = sp.Rational(8*14*11, 17**3)
    assert norm_obstruction_at(shortest, 3)
    assert norm_obstruction_at(complement, 7)
    triple = [sp.Matrix(v) for v in ((0, 1, 4), (0, 4, 1), (2, 2, 3))]
    assert [p.dot(q) for p, q in combinations(triple, 2)] == [8, 14, 11]

    # An odd permutation of spatial coordinates reverses the signed loop.
    # On two identical copies this map commutes with the independent clock.
    c = multiply_by_vector(sp.Matrix([1, -4, 1]), identity, right=True)
    for i, j in combinations(range(3), 2):
        n = sp.eye(3)[:, i] - sp.eye(3)[:, j]
        swap = sp.eye(3) - n*n.T
        lift = multiply_by_vector(n, identity, right=False)*c/6
        assert lift.T * lift == identity and lift*right*lift.T == -right
        for index, p in enumerate(points):
            assert lift*planes[index]*lift.T == planes[points.index(swap*p)]

    ijk = [points.index(sp.Matrix(v)) for v in ((3, 2, 2), (2, 3, 2), (2, 2, 3))]
    p, q, r = [planes[t] for t in ijk]
    loop = sp.diag(p*q*r*p, p*q*r*p)
    clock = sp.Matrix([[0, -1], [1, 0]])
    clock = sp.kronecker_product(clock, identity)
    handed = sp.diag(right, right)*clock
    assert handed.T == handed and handed**2 == 17*sp.eye(8)
    for lam in (sp.Rational(1, 5), sp.Rational(-1, 5)):
        density = sp.eye(8) + lam*handed
        assert density.is_positive_definite
        assert (density*loop).trace()/4 == sp.Rational(65, 68)
        assert (density*loop*clock).trace()/4 == 7*lam/68

    data = json.loads(Path(__file__).with_name("ch39-five-tick-record.json").read_text())
    den = 2**data["bits"]
    assert len(data["rows"]) == 46 and sum(row[2] for row in data["rows"]) == den
    for offset, target in enumerate(data["levels"]):
        total = sine = 0
        for _, height, weight, (sign, chart), _ in data["rows"]:
            assert weight > 0 and height > 0 and sign in (-1, 1)
            mass = Fraction(weight * height**2, height**2 + (offset*den)**2)
            turn = Fraction(2*sign*chart*den, den**2 + chart**2)
            total += mass // 1
            sine += (mass*turn) // 1
        assert max(abs(total-target-1), abs(total+46-target)) < 128
        if offset == 0:
            factor = Fraction(14, 340*den)
            assert Fraction("-0.0170319270") < sine*factor <= (sine+46)*factor < Fraction("-0.0170319269")
    print("81 compressions, 504 oriented triples, 3 coordinate transpositions and 42 arithmetic offsets verified.")
    print("At the normalization offset the signed difference lies in (-0.0170319270, -0.0170319269).")


if __name__ == "__main__":
    main()

