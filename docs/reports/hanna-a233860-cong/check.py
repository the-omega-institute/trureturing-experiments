#!/usr/bin/env python3
"""Exact finite checks for all three COMMENTS conjectures of OEIS A233860.

Only Python's integer arithmetic and fractions.Fraction are used as evidence.
These finite checks do not prove the unbounded conjectures.
"""
from fractions import Fraction
from math import factorial

N = 160
M = 32
SOURCE_FIRST_20 = [
    1, 1, 5, 35, 355, 4465, 69125, 1252475, 26151475, 616872025,
    16234589525, 471382586675, 14970245087875, 516142537458625,
    19199498482905125, 766394702651760875, 32676482018991377875,
    1482055899582130035625, 71248344993651091083125,
    3618867148116847594611875,
]
EXCEPTION_VALUES = {
    3: 5, 4: 35, 8: 1252475, 9: 26151475, 13: 14970245087875,
    14: 516142537458625, 18: 1482055899582130035625,
    19: 71248344993651091083125, 23: 641829656892378802582108653125,
    24: 39577486974389791588982993046875,
}


def mul(a, b, degree, modulus=None):
    c = [0] * (degree + 1)
    for i, ai in enumerate(a[:degree + 1]):
        if ai:
            for j, bj in enumerate(b[:degree + 1 - i]):
                if bj:
                    c[i + j] += ai * bj
    return c if modulus is None else [x % modulus for x in c]


def compose(f, g, degree, modulus=None):
    result = [0] * (degree + 1)
    for c in reversed(f[:degree + 1]):
        result = mul(result, g, degree, modulus)
        result[0] += c
        if modulus is not None:
            result[0] %= modulus
    return result


def valuation(n, p):
    assert n != 0
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count


def main():
    # P: independent unrestricted-partition dynamic program.
    partition = [1] + [0] * N
    for part in range(1, N + 1):
        for n in range(part, N + 1):
            partition[n] += partition[n - part]

    # E: finite Euler product; factors above N cannot affect this truncation.
    euler = [1] + [0] * N
    for part in range(1, N + 1):
        for n in range(N, part - 1, -1):
            euler[n] -= euler[n - part]
    assert mul(euler, partition, N) == [1] + [0] * N

    # Compare the product with generalized pentagonal exponents, both signs.
    pentagonal = [0] * (N + 1)
    for j in range(-N, N + 1):
        r = j * (3 * j - 1) // 2
        if r <= N:
            pentagonal[r] += 1 if j % 2 == 0 else -1
    assert euler == pentagonal

    # Integer derivation D = P d/dt. Each iteration loses one known degree.
    # a(n) = (D^(n-1)P)(0), calculated without rational series reversion.
    current = partition[:N]
    a = [0]
    for n in range(1, N + 1):
        a.append(current[0])
        degree = len(current) - 2
        if degree >= 0:
            deriv = [(i + 1) * current[i + 1] for i in range(degree + 1)]
            current = mul(partition, deriv, degree)
    assert a[1:21] == SOURCE_FIRST_20
    print("FIDELITY: OEIS A233860 revision 10 DATA, all 20 terms reproduced")
    print("FIRST_20 =", a[1:21])
    print("EULER_PRODUCT: reciprocal and pentagonal coefficients agree through degree", N)

    # G = integral E, with zero integration constant.
    primitive = [Fraction(0)] + [Fraction(euler[r], r + 1) for r in range(N)]
    # Recursive reversion, independent of the integer derivation calculation.
    inverse = [Fraction(0)] * (M + 1)
    inverse[1] = Fraction(1)
    for n in range(2, M + 1):
        inverse[n] = -compose(primitive, inverse, n)[n]
    assert compose(primitive, inverse, M) == [0, 1] + [0] * (M - 1)
    assert compose(inverse, primitive, M) == [0, 1] + [0] * (M - 1)
    assert [factorial(n) * inverse[n] for n in range(M + 1)] == a[:M + 1]
    ode_rhs = compose(partition, inverse, M - 1)
    assert ode_rhs == [(n + 1) * inverse[n + 1] for n in range(M)]
    print("INDEPENDENT_RATIONAL_INVERSE: both compositions through degree", M,
          "; source ODE through degree", M - 1, "; factorial-scaled values agree")

    for p in (5, 7):
        residues = {j * (3 * j - 1) // 2 % p for j in range(2 * p)}
        assert p - 1 not in residues
        assert all(c.denominator % p != 0 for c in primitive)
        assert all(c.denominator % p != 0 for c in inverse)
        assert all(a[n] % p ** valuation(factorial(n), p) == 0
                   for n in range(1, N + 1))
        print("PRIME_LOCAL:", p, "pentagonal residues", sorted(residues),
              "; G denominators checked through", N,
              "; inverse denominators through", M,
              "; factorial lower bound through", N)

    g_sparse = {1: 1, 2: 2, 3: 3, 6: 1, 8: 2, 13: 3, 16: 4, 23: 2}
    b_sparse = {1: 1, 2: 3, 5: 4, 6: 2, 10: 3, 11: 3, 12: 3,
                15: 4, 16: 3, 17: 1, 20: 1, 21: 4, 22: 2}
    g = [g_sparse.get(i, 0) for i in range(25)]
    b = [b_sparse.get(i, 0) for i in range(25)]
    reduce5 = lambda q: q.numerator * pow(q.denominator, -1, 5) % 5
    assert g == [reduce5(q) for q in primitive[:25]]
    assert b == [reduce5(q) for q in inverse[:25]]
    certificate = compose(g, b, 24, 5)
    assert certificate == [0, 1] + [0] * 23
    assert compose(b, g, 24, 5) == certificate
    print("F5_G_SPARSE =", g_sparse)
    print("F5_B_SPARSE =", b_sparse)
    print("F5_G_OF_B_DEGREES_0_TO_24 =", certificate)

    exceptional = [n for n in range(1, 25)
                   if valuation(factorial(n), 5) < (n + 2) // 5]
    assert exceptional == sorted(EXCEPTION_VALUES)
    for n in exceptional:
        assert b[n] == 0
        assert a[n] == EXCEPTION_VALUES[n]
        assert valuation(a[n], 5) >= valuation(factorial(n), 5) + 1
        print("EXCEPTION:", n, "a(n) =", a[n], "v5 =", valuation(a[n], 5))
    assert all(valuation(factorial(n), 5) >= (n + 2) // 5
               for n in range(25, N + 1))

    five_cases = seven_cases = 0
    for n in range(1, N + 1):
        assert a[n] % 2 == 1
        assert valuation(a[n], 5) >= (n + 2) // 5
        assert valuation(a[n], 7) >= n // 7
        for k in range(1, (n + 2) // 5 + 1):
            assert n >= 5 * k - 2 and a[n] % 5 ** k == 0
            five_cases += 1
        for k in range(1, n // 7 + 1):
            assert n >= 7 * k and a[n] % 7 ** k == 0
            seven_cases += 1
    print("CONJECTURES: all n = 1..", N, "; parity cases =", N,
          "; admissible (k,n) for 5 =", five_cases,
          "; admissible (k,n) for 7 =", seven_cases)
    print("SCOPE: finite exact checks and the displayed polynomial certificate;")
    print("       no unbounded theorem or Lean proof is certified by this program.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
