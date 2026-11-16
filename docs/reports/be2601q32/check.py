#!/usr/bin/env python3
"""Exact bounded checks for arXiv:2601.07057v3, Theorem 3.1 / Question 3.2.

The finite sets below are support windows in infinite groups, not quotient
quandles. Only coefficients are reduced modulo primes. No general theorem
is certified by this program. Standard-library code, independently written.
"""
from fractions import Fraction
from itertools import combinations, product


def add_core(x, y):
    if isinstance(x, tuple):
        return tuple(2 * b - a for a, b in zip(x, y))
    return 2 * y - x


def hmul(x, y):
    a, b, c = x
    d, e, f = y
    return a + d, b + e, c + f + a * e


def hinv(x):
    a, b, c = x
    return -a, -b, -c + a * b


def hcore(x, y):
    return hmul(hmul(y, hinv(x)), y)


def square(u, core, modulus=None):
    answer = {}
    for x, a in u.items():
        for y, b in u.items():
            z = core(x, y)
            answer[z] = answer.get(z, 0) + a * b
            if modulus is not None:
                answer[z] %= modulus
    return {x: a for x, a in answer.items() if a != 0}


def verify(u, core, modulus=None):
    sq = square(u, core, modulus)
    trivial = not u or (len(u) == 1 and next(iter(u.values())) == 1)
    assert (sq == u) == trivial, (u, sq)
    if len(u) >= 2:
        a, b = min(u), max(u)
        high, low = core(a, b), core(b, a)
        assert high > b and low < a, (a, b, high, low)
        high_pairs = [(x, y) for x in u for y in u if core(x, y) == high]
        low_pairs = [(x, y) for x in u for y in u if core(x, y) == low]
        assert high_pairs == [(a, b)]
        assert low_pairs == [(b, a)]
        assert all(low <= core(x, y) <= high for x in u for y in u)
        high_coefficient = u[a] * u[b]
        low_coefficient = u[b] * u[a]
        if modulus is not None:
            high_coefficient %= modulus
            low_coefficient %= modulus
        assert high_coefficient != 0 and low_coefficient != 0
        assert sq[high] == high_coefficient and sq[low] == low_coefficient
        assert high not in u and low not in u
    return int(sq == u)


def full_window(label, points, coefficients, modulus=None):
    cases = idempotents = 0
    for values in product(coefficients, repeat=len(points)):
        u = {x: a for x, a in zip(points, values) if a != 0}
        idempotents += verify(u, add_core, modulus)
        cases += 1
    assert idempotents == len(points) + 1
    print(f"{label}: vectors={cases}; idempotents={idempotents}; "
          "classification and both unique extremes OK")


def sparse_window(label, points, coefficients, core, modulus=None):
    cases = idempotents = 0
    for size in range(4):
        for support in combinations(points, size):
            for values in product(coefficients, repeat=size):
                u = dict(zip(support, values))
                idempotents += verify(u, core, modulus)
                cases += 1
    assert idempotents == len(points) + 1
    print(f"{label}: supports<=3; vectors={cases}; idempotents={idempotents}; "
          "classification and both unique extremes OK")


def heisenberg_structure(points):
    identity = (0, 0, 0)
    for x in points:
        assert hmul(x, hinv(x)) == identity == hmul(hinv(x), x)
        assert hcore(x, x) == x
        for y in points:
            assert hcore(hcore(x, y), y) == x
            for z in points:
                assert hmul(hmul(x, y), z) == hmul(x, hmul(y, z))
                assert hcore(hcore(x, y), z) == hcore(hcore(x, z), hcore(y, z))
                if x < y:
                    assert hmul(z, x) < hmul(z, y)
                    assert hmul(x, z) < hmul(y, z)
    x, y = (1, 0, 0), (0, 1, 0)
    assert hmul(x, y) != hmul(y, x)
    print(f"Heisenberg structure: {len(points) ** 3} triples; associativity, "
          "inverse, both order translations and right core laws OK; noncommutative")


def additive_identity(points):
    cases = 0
    for a, b in combinations(points, 2):
        for x in points:
            for y in points:
                if a <= x <= b and a <= y <= b:
                    gap = (2 * b - a) - (2 * y - x)
                    assert gap == 2 * (b - y) + (x - a)
                    assert gap >= 0
                    assert (gap == 0) == (x == a and y == b)
                    cases += 1
    print(f"Additive gap identity: {cases} in-range (a,b,x,y) tuples OK")


def main():
    # Fidelity: exactly the product and extreme formula in Theorem 3.1.
    integers = list(range(-2, 3))
    full_window("Fidelity, G=Z, k=Z, coefficients=-2..2", integers, integers)
    additive_identity(list(range(-4, 5)))
    # The expected extension, with index doubling always in the infinite group.
    rational_coefficients = [Fraction(-1), Fraction(-1, 2), Fraction(0),
                             Fraction(1, 2), Fraction(1)]
    full_window("G=Z, k=Q", integers, rational_coefficients)
    for p in (2, 3, 5):
        full_window(f"G=Z, k=F_{p}", integers, list(range(p)), p)
    rank_two = list(product(range(-1, 2), repeat=2))
    sparse_window("G=Z^2 lex, k=Z", rank_two, [-2, -1, 1, 2], add_core)
    sparse_window("G=Z^2 lex, k=Q", rank_two,
                  [Fraction(-1), Fraction(-1, 2), Fraction(1, 2), Fraction(1)], add_core)
    for p in (2, 3, 5):
        sparse_window(f"G=Z^2 lex, k=F_{p}", rank_two, list(range(1, p)), add_core, p)
    heisenberg = list(product(range(-1, 2), repeat=3))
    heisenberg_structure(heisenberg)
    sparse_window("G=Heisenberg(Z) lex, k=Z", heisenberg, [-1, 1], hcore)
    sparse_window("G=Heisenberg(Z) lex, k=F_2", heisenberg, [1], hcore, 2)
    sparse_window("G=Heisenberg(Z) lex, k=F_3", heisenberg, [1, 2], hcore, 3)
    # The singleton coefficient step in the argument is also checked separately.
    rationals = sorted({Fraction(a, b) for a in range(-4, 5) for b in range(1, 5)})
    for a in rationals:
        assert (a * a == a) == (a in (0, 1))
    print(f"Singleton coefficient equation: {len(rationals)} exact rational scalars OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
