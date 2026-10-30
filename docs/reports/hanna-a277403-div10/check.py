#!/usr/bin/env python3
"""Exact finite checks for Paul D. Hanna's OEIS A277403 conjecture.

Source: https://oeis.org/A277403, revision #21, offset 1.
Uses only Python's standard library; no floating point or external packages.
Finite checks support a preregistration, not an all-index formal proof.
"""
from fractions import Fraction
from math import comb, factorial

SOURCE_TERMS = (
    1, 2, 10, 90, 1190, 20930, 462070, 12326790, 386855630,
    14000898310, 575440398330, 26532920708070, 1358954912773010,
    76682330257445570, 4734315243483414890, 317932511564758225170,
    23106045191162625194230, 1809303767549542227341490,
    152057767850058496005946030, 13668688227104664304597942910,
    1310201986290043690952261887230, 133552478071366935949713096470670,
    14440878313638992240490923468851610,
)


def bell_solution(bound, modulus=None):
    """Solve coefficient equations using partial exponential Bell polynomials."""
    a = [0] * (bound + 1)
    bell = [[0] * (bound + 1) for _ in range(bound + 1)]
    bell[0][0] = 1
    for n in range(1, bound + 1):
        # x - integral(A): b_1 = 1, b_j = -a_{j-1} for j >= 2.
        for k in range(1, n + 1):
            value = sum(
                comb(n - 1, j - 1)
                * (1 if j == 1 else -a[j - 1])
                * bell[n - j][k - 1]
                for j in range(1, n - k + 2)
            )
            bell[n][k] = value if modulus is None else value % modulus
        assert bell[n][n] == 1
        value = (n == 1) + a[n - 1] - sum(
            a[k] * bell[n][k] for k in range(1, n)
        )
        a[n] = value if modulus is None else value % modulus
        lhs = sum(a[k] * bell[n][k] for k in range(1, n + 1))
        rhs = (n == 1) + a[n - 1]
        assert lhs == rhs if modulus is None else (lhs - rhs) % modulus == 0
    return a


def mul(f, g, bound):
    """Ordinary polynomial multiplication over Fraction, truncated in degree."""
    out = [Fraction(0) for _ in range(bound + 1)]
    for i, x in enumerate(f[: bound + 1]):
        if x:
            for j, y in enumerate(g[: bound + 1 - i]):
                if y:
                    out[i + j] += x * y
    return out


def compose(f, g, bound):
    assert g[0] == 0
    out = [Fraction(0) for _ in range(bound + 1)]
    for c in reversed(f[: bound + 1]):
        out = mul(out, g, bound)
        out[0] += c
    return out


def integral(f, bound):
    return [Fraction(0)] + [Fraction(f[n - 1], n) for n in range(1, bound + 1)]


def rational_solution(bound):
    """Independent ordinary-power-series substitution, without Bell polynomials."""
    a = [Fraction(0) for _ in range(bound + 1)]
    for n in range(1, bound + 1):
        g = integral(a, bound)
        b = [-c for c in g]
        b[1] += 1
        c = compose(a, b, n)
        # The not-yet-set coefficient a[n] has multiplier b[1]**n = 1.
        a[n] = Fraction(n == 1) + g[n] - c[n]
    g = integral(a, bound)
    b = [-c for c in g]
    b[1] += 1
    rhs = g.copy()
    rhs[1] += 1
    assert compose(a, b, bound) == rhs
    return [factorial(n) * c for n, c in enumerate(a)]


def main():
    bound = 60
    a = bell_solution(bound)
    assert tuple(a[1 : len(SOURCE_TERMS) + 1]) == SOURCE_TERMS
    print('FIDELITY: all 23 published %S/%T/%U terms reproduced exactly (n=1..23)')
    print('FIRST_TEN:', ','.join(map(str, a[1:11])))
    independent = rational_solution(24)
    assert independent == a[:25]
    assert all(c.denominator == 1 for c in independent)
    print('INDEPENDENT: ordinary Fraction substitution agrees with integer Bell recursion for n=0..24')
    failures = [n for n in range(3, bound + 1) if a[n] % 10]
    assert failures == []
    print('DIVISIBILITY: n=3..60, failures=[]')
    for modulus in (2, 5, 10):
        residues = bell_solution(bound, modulus)
        assert residues == [c % modulus for c in a]
        assert residues[1:] == [1 % modulus, 2 % modulus] + [0] * (bound - 2)
        print(f'MODULAR: recursion modulo {modulus} agrees with integer reduction and polynomial witness through n=60')
    # Compute the entire polynomial residual over Q, not in Z/10Z.
    degree = 6
    p = list(map(Fraction, [0, 1, 1, 0, 0, 0, 0]))
    j = integral(p, degree)
    b = [-c for c in j]
    b[1] += 1
    residual = compose(p, b, degree)
    residual = [residual[n] - j[n] - (n == 1) for n in range(degree + 1)]
    expected = [Fraction(0), Fraction(0), Fraction(0), Fraction(-5, 3),
                Fraction(-5, 12), Fraction(1, 3), Fraction(1, 9)]
    assert residual == expected
    ecoeff = [factorial(n) * c for n, c in enumerate(residual)]
    assert ecoeff == [0, 0, 0, -10, -10, 40, 80]
    assert all(c.denominator == 1 and c.numerator % 10 == 0 for c in ecoeff)
    # deg(P)=2 and deg(x-J)=3, hence no residual terms above degree 6.
    print('EXACT_RESIDUAL_ORDINARY:', ','.join(map(str, residual)))
    print('EXACT_RESIDUAL_EGF:', ','.join(str(int(c)) for c in ecoeff))
    print('DEGREE: deg(P)=2, deg(x-J)=3, so residual coefficients above degree 6 are zero')
    print('SCOPE: finite arithmetic checks and exact polynomial identity; no Lean proof of the conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
