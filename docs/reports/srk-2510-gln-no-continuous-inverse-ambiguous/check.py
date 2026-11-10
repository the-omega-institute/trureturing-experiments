#!/usr/bin/env python3
"""Exact algebra checks; no finite test certifies the missing topology of GL_n(C)."""

from fractions import Fraction
from itertools import permutations, product


def fidelity_finite_groups():
    """Source Theorem 1.3: cyclic groups of orders 1 through 8, all bijections."""
    rows = []
    for m in range(1, 9):
        count = 0
        for p in permutations(range(m)):
            square = all(p[p[x]] == (-x) % m for x in range(m))
            pinv = [0] * m
            for x, y in enumerate(p):
                pinv[y] = x
            ambiguous = all(pinv[x] == (-p[x]) % m for x in range(m))
            assert square == ambiguous
            count += square
        nonselfinverse = sum((2 * x) % m != 0 for x in range(m))
        assert (count != 0) == (nonselfinverse % 4 == 0)
        rows.append((m, nonselfinverse, count))
    print("FIDELITY Theorem 1.3: (order, nonselfinverse, number_of_roots)", rows)


def fidelity_even_tori():
    """Source section 2 construction, restricted to roots of unity as exact exponents."""
    cases = 0
    for m in range(2, 10):
        for dimension in (2, 4):
            points = list(product(range(m), repeat=dimension))

            def f(x):
                return tuple(y for i in range(0, dimension, 2)
                             for y in (x[i + 1], (-x[i]) % m))

            assert len({f(x) for x in points}) == len(points)
            for x in points:
                assert f(f(x)) == tuple((-z) % m for z in x)
                cases += 1
    print("FIDELITY even-torus formula: exact torsion points checked", cases)


def det(a):
    b = [[Fraction(x) for x in row] for row in a]
    n = len(b)
    ans = Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if b[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            b[j], b[pivot] = b[pivot], b[j]
            ans = -ans
        v = b[j][j]
        ans *= v
        for i in range(j + 1, n):
            ratio = b[i][j] / v
            for k in range(j, n):
                b[i][k] -= ratio * b[j][k]
    return ans


def inverse(a):
    n = len(a)
    b = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if b[i][j])
        b[j], b[pivot] = b[pivot], b[j]
        v = b[j][j]
        b[j] = [x / v for x in b[j]]
        for i in range(n):
            if i != j:
                v = b[i][j]
                b[i] = [x - v * y for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b]


def matrix_checks():
    """Rational matrices embed in complex matrices; check inversion/determinant identities."""
    cases = 0
    for n in range(2, 9):
        for t in range(-5, 6):
            if t == 0:
                continue
            a = [[Fraction(i == j) for j in range(n)] for i in range(n)]
            a[0][0] = Fraction(t)
            # Row additions preserve determinant and ensure off-diagonal examples occur.
            for i in range(1, n):
                a[i] = [x + (i + t) * y for x, y in zip(a[i], a[i - 1])]
            ai = inverse(a)
            assert det(a) == t
            assert det(ai) == Fraction(1, t)
            assert inverse(ai) == a
            assert [[sum(a[i][k] * ai[k][j] for k in range(n))
                     for j in range(n)] for i in range(n)] == [
                         [Fraction(i == j) for j in range(n)] for i in range(n)]
            cases += 1
    print("SETTLEMENT ingredient: det(A^-1)=det(A)^-1 on exact GL_n(Q) examples", cases)


def rank_one_obstruction():
    # Every integer has one of these residues. This is an exhaustive congruence
    # certificate, rather than evidence from a bounded search for k.
    squares = {r * r % 4 for r in range(4)}
    assert squares == {0, 1}
    assert (-1) % 4 not in squares
    for k in range(-100, 101):
        endomorphism = lambda x: k * x
        for x in range(-10, 11):
            assert endomorphism(endomorphism(x)) == (k * k) * x
        assert endomorphism(endomorphism(1)) != -1
    print("SETTLEMENT rank-one obstruction: squares modulo 4 are", sorted(squares))


def insufficient_surjection_check():
    # A split surjection to Z is insufficient without invariance of its kernel.
    # On Z^2, J^2=-I, but the first projection does not descend through J.
    j = lambda p: (p[1], -p[0])
    for x, y in product(range(-5, 6), repeat=2):
        assert j(j((x, y))) == (-x, -y)
    assert (0, 0)[0] == (0, 1)[0]
    assert j((0, 0))[0] != j((0, 1))[0]
    print("ROUTE CHECK: a split surjection Z^2 -> Z alone does not give rank-one functoriality")


def main():
    fidelity_finite_groups()
    fidelity_even_tori()
    matrix_checks()
    rank_one_obstruction()
    insufficient_surjection_check()
    print("LIMIT: topology and naturality for GL_n(C) remain unverified; this is not a settlement")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
