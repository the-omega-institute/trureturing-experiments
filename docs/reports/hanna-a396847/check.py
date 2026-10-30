#!/usr/bin/env python3
"""Exact finite checks for Paul D. Hanna's OEIS A396847, revision 18.

Only integers and fractions.Fraction carry evidence. These checks support a
preregistration; they do not prove the infinite boundary-convergence claim.
Source: https://oeis.org/A396847/internal (author July 12, 2026).
"""
from fractions import Fraction
from itertools import combinations


def primes_through(bound):
    result = []
    for n in range(2, bound + 1):
        if all(n % p for p in result if p * p <= n):
            result.append(n)
    return result


def mul(a, b, degree):
    c = [0] * (degree + 1)
    for i, x in enumerate(a[:degree + 1]):
        if x:
            for j, y in enumerate(b[:degree + 1 - i]):
                c[i + j] += x * y
    return c


def power(a, exponent, degree):
    result = [1] + [0] * degree
    base = a
    while exponent:
        if exponent % 2:
            result = mul(result, base, degree)
        exponent //= 2
        if exponent:
            base = mul(base, base, degree)
    return result


def step(a, exponent, degree):
    increment = power(a, exponent, degree)
    return [x + y for x, y in zip(a, increment)]


def evaluate(a, x):
    result = Fraction(0)
    for coefficient in reversed(a):
        result = result * x + coefficient
    return result


def rounded_orbit(numerator, exponents, upward):
    scale = 10 ** 30
    m = numerator * scale // 1000
    for p in exponents:
        denominator = scale ** (p - 1)
        term, remainder = divmod(m ** p, denominator)
        if upward and remainder:
            term += 1
        m += term
    return m


def main():
    primes = primes_through(61)
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19, 23,
                      29, 31, 37, 41, 43, 47, 53, 59, 61]
    degree = 40
    a = [0, 1] + [0] * (degree - 1)
    snapshots = []
    for p in primes:
        old = a
        a = step(a, p, degree)
        assert a[0] == 0 and a[1] == 1
        assert all(x <= y for x, y in zip(old, a))
        assert a[:min(p, degree + 1)] == old[:min(p, degree + 1)]
        snapshots.append(a)
    published = [1, 1, 1, 3, 4, 6, 16, 52, 148, 364, 861, 1983,
                 4458, 10260, 24881, 63339, 167035, 448673, 1207450,
                 3226010, 8544369, 22478369, 58905246, 154187918,
                 403928406, 1060043106, 2787736717, 7348365111,
                 19422984322, 51509520556, 137150272684, 366847701440,
                 986008997934, 2662974001118, 7224478109771, 19678015523249]
    assert a[1:len(published) + 1] == published
    assert snapshots[1][:7] == [0, 1, 1, 1, 3, 3, 1]
    print('FIDELITY: all 36 published coefficients and the full F_2 expansion agree')
    # After the exponent exceeds K, coefficients of degrees <= K cannot change.
    for k in range(degree + 1):
        stable_stage = sum(p <= k for p in primes)
        base = ([0, 1] + [0] * (degree - 1)) if stable_stage == 0 else snapshots[stable_stage - 1]
        assert all(s[:k + 1] == base[:k + 1] for s in snapshots[stable_stage:])
    print('STABILIZATION: degrees 0..40 checked through exponent 61')
    # Independently compare full polynomial evaluation and scalar recursion.
    a = [0, 1] + [0] * 29
    for n, p in enumerate(primes[:3], 1):
        a = step(a, p, 30)
        for x in [Fraction(0), Fraction(1, 6), Fraction(1, 3), Fraction(1, 2)]:
            y = x
            for q in primes[:n]:
                y += y ** q
            assert evaluate(a, x) == y
    print('BRIDGE: full F_1, F_2, F_3 evaluation equals scalar recursion at four rationals')
    cases = 0
    for length in range(1, 4):
        for exponents in combinations(range(2, 9), length):
            y = Fraction(1, 6)
            for p in exponents:
                y += y ** p
                assert y <= Fraction(1, 3)
            cases += 1
    print('SMALL CASES: lower trapping invariant holds for', cases, 'increasing exponent lists')
    scale = 10 ** 30
    upper = rounded_orbit(356, primes[:15], upward=True)
    lower = rounded_orbit(357, primes[:11], upward=False)
    assert upper == 861701822816469940744335134135
    assert lower == 1653493773472349502191991892877
    assert Fraction(upper, scale) < Fraction(862, 1000)
    assert lower > scale
    b = Fraction(87, 100)
    tail = b ** 53 / (1 - b)
    assert 200 * 87 ** 53 < 13 * 100 ** 52
    assert tail < Fraction(1, 200)
    assert Fraction(862, 1000) + Fraction(1, 200) < b
    print('UPPER ORBIT: F_15(356/1000) <=', str(upper) + '/' + str(scale), '< 862/1000')
    print('TAIL: (87/100)^53 / (1 - 87/100) < 1/200; 862/1000 + 1/200 < 87/100')
    print('LOWER ORBIT: F_11(357/1000) >=', str(lower) + '/' + str(scale), '> 1')
    print('RADIUS BRACKET: 356/1000 < R < 357/1000, conditional on the stated series/orbit bridge')
    print('No assertion certifies the longer decimal printed in OEIS or proves the infinite claim.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
