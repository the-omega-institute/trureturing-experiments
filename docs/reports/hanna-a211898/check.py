#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A211898 valuation conjecture.

Primary source: https://oeis.org/A211898/internal, revision 9.
Finite checks support a proposed universal argument; they are not its proof.
"""
from fractions import Fraction
from math import comb, factorial

SOURCE_PREFIX = (
    1, 3, 9, 261, 13419, 7867287, 10444212819, 84955235950827,
    2235017786095822257, 273416315791427558035965,
    125533366255776787874473759857,
    242979442003484538229530424638338553,
    1852958949086213206247388599213928431454549,
)
N = 120
DIRECT_N = 24


def v3(x):
    """Ordinary finite 3-adic valuation of a positive integer."""
    assert isinstance(x, int) and x > 0
    result = 0
    while x % 3 == 0:
        x //= 3
        result += 1
    return result


def mul_truncated(f, g, degree):
    return [sum((f[k] * g[n-k] for k in range(n+1)), Fraction(0))
            for n in range(degree+1)]


def main():
    d = [2**k - (-1)**k for k in range(N+1)]
    b = [d[k]**k for k in range(N+1)]
    # Original logarithmic recurrence, obtained from A' = A H'.
    a = [Fraction(1)]
    for n in range(1, N+1):
        a.append(sum((b[k] * a[n-k] for k in range(1, n+1)),
                     Fraction(0)) / n)
    assert tuple(a[:len(SOURCE_PREFIX)]) == SOURCE_PREFIX

    # Independent coefficient extraction from exp(H) = sum H^j/j!.
    h = [Fraction(0)] + [Fraction(b[k], k) for k in range(1, DIRECT_N+1)]
    h_power = [Fraction(1)] + [Fraction(0)] * DIRECT_N
    direct = [Fraction(0)] * (DIRECT_N+1)
    for j in range(DIRECT_N+1):
        for n in range(DIRECT_N+1):
            direct[n] += h_power[n] / factorial(j)
        h_power = mul_truncated(h_power, h, DIRECT_N)
    assert direct == a[:DIRECT_N+1]

    assert all(d[k] % 3 == 0 for k in range(1, N+1))
    c = [x // 3 for x in d]
    assert c[1:4] == [1, 1, 3]
    # Independent integer, factorial-scaled recurrence.
    r = [1]
    for n in range(1, N+1):
        r.append(sum(comb(n-1, k-1) * factorial(k-1) * c[k]**k * r[n-k]
                     for k in range(1, n+1)))
    involutions = [1, 1]
    for n in range(2, N+1):
        involutions.append(involutions[-1] + (n-1)*involutions[-2])

    for n in range(N+1):
        assert a[n].denominator == 1 and a[n] > 0
        assert Fraction(factorial(n), 3**n) * a[n] == r[n]
        assert r[n] % 3 == (1, 1, 2)[n % 3]
        assert r[n] % 3 == involutions[n] % 3
        assert r[n] % 3 != 0
        assert v3(a[n].numerator) == n - v3(factorial(n))
        if n > 0:
            assert a[n].numerator % 6 == 3
        if n >= 2:
            assert (r[n] - r[n-1] - (n-1)*r[n-2]) % 3 == 0
    for k in range(3, N+1):
        assert factorial(k-1) * c[k]**k % 3 == 0

    print(f"FIDELITY: all {len(SOURCE_PREFIX)} source coefficients agree (indices 0..12)")
    print(f"DIRECT_EXP: exact truncated exponential agrees (indices 0..{DIRECT_N})")
    print(f"SCALED: independent recurrence and n!*a(n)/3^n agree (indices 0..{N})")
    print(f"VALUATION: positive integral coefficients and v3(a(n)) = n-v3(n!) (indices 0..{N})")
    print(f"RESIDUES: scaled coefficients and involution numbers repeat 1,1,2 (indices 0..{N})")
    print(f"SOURCE_FORMULA: a(n) = 3 mod 6 (indices 1..{N})")
    print("a(0..5) =", ", ".join(str(x.numerator) for x in a[:6]))
    print("r(0..5) =", ", ".join(str(x) for x in r[:6]))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
