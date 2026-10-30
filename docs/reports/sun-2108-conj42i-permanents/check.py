#!/usr/bin/env python3
"""Exact finite checks for Sun, arXiv:2108.07723v7, Conjecture 4.2(i).

Fidelity: Remark 4.2's two q-integer determinant identities.
Target: both absolute-value permanent congruences, with independent algorithms.
All arithmetic used as evidence is integer arithmetic. No external packages.
"""
from itertools import permutations, product
from math import prod


def prime(n):
    return n >= 2 and all(n % d for d in range(2, n) if d * d <= n)


def matrix(n, shift):
    # The source uses labels 1,...,n, and takes absolute values in Z.
    return [[abs(j - k + shift) for k in range(1, n + 1)]
            for j in range(1, n + 1)]


def permanent_dp(a):
    n = len(a)
    dp = [0] * (1 << n)
    dp[0] = 1
    for mask in range(1, 1 << n):
        row = a[mask.bit_count() - 1]
        bits = mask
        total = 0
        while bits:
            bit = bits & -bits
            total += row[bit.bit_length() - 1] * dp[mask ^ bit]
            bits ^= bit
        dp[mask] = total
    return dp[-1]


def permanent_brute(a):
    return sum(prod(a[i][s[i]] for i in range(len(a)))
               for s in permutations(range(len(a))))


def permanent_ryser(a):
    # Inclusion-exclusion is independent of the row-by-row subset recurrence.
    n = len(a)
    total = 0
    for mask in range(1 << n):
        cols = [j for j in range(n) if mask & (1 << j)]
        total += (-1) ** (n - len(cols)) * prod(
            sum(row[j] for j in cols) for row in a)
    return total


def determinant_bareiss(a):
    a = [row[:] for row in a]
    n = len(a)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((r for r in range(k, n) if a[r][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def fidelity():
    count = 0
    for n in range(2, 10):
        for q in (0, 1, 2, 3):
            def qinteger(v):
                return sum(q ** t for t in range(v))
            for shift in (0, 1):
                a = [[qinteger(v) for v in row] for row in matrix(n, shift)]
                expected = ((-1) ** (n - 1) * (n - 1) if shift == 0 else 1)
                expected *= (1 + q) ** (n - 2)
                actual = determinant_bareiss(a)
                assert actual == expected, (n, q, shift, actual, expected)
                count += 1
    print(f'FIDELITY Remark 4.2: {count} exact determinant identities; '
          'n=2..9, q=0,1,2,3, shifts=0,1')
    # An entirely separate determinant algorithm cross-checks the implementation.
    for n in range(2, 7):
        for shift in (0, 1):
            a = matrix(n, shift)
            value = sum((-1) ** sum(s[i] > s[j]
                                    for i in range(n) for j in range(i + 1, n))
                        * prod(a[i][s[i]] for i in range(n))
                        for s in permutations(range(n)))
            assert value == determinant_bareiss(a)
    print('FIDELITY_CROSSCHECK Leibniz determinants agree for n=2..6')


def check_targets():
    for p in (p for p in range(3, 20) if prime(p)):
        values = []
        for shift in (0, 1):
            a = matrix(p, shift)
            value = permanent_dp(a)
            assert (2 * value + (1 if shift == 0 else -1)) % p == 0
            if p <= 7:
                assert value == permanent_brute(a)
            if p <= 13:
                assert value == permanent_ryser(a)
            values.append(value)
        print(f'TARGET p={p}: per(A0)={values[0]}; per(A1)={values[1]}; '
              f'residues={values[0] % p},{values[1] % p}; both congruences OK')
    print('PERMANENT_CROSSCHECK direct permutation sums p<=7; Ryser p<=13')


def check_reduction():
    count = 0
    for p in (p for p in range(3, 102) if prime(p)):
        for shift in (0, 1):
            for i in range(p):
                for j in range(p):
                    r = (i - j) % p
                    lhs = (-1) ** (i + j) * abs(i - j + shift)
                    rhs = (-1) ** r * (r + shift)
                    assert (lhs - rhs) % p == 0, (p, shift, i, j)
                    count += 1
            total = sum((-1) ** r * (r + shift) for r in range(p))
            assert total == (p - 1) // 2 + shift
    print(f'SIGN_REDUCTION {count} entries; all odd primes p<=101, shifts=0,1')


def check_classical_ingredient():
    count = 0
    # Exhaust every function F_p -> F_p for p=3,5.
    for p in (3, 5):
        for f in product(range(p), repeat=p):
            a = [[f[(i - j) % p] for j in range(p)] for i in range(p)]
            assert (permanent_dp(a) - sum(f)) % p == 0
            count += 1
    print(f'CIRCULANT_INGREDIENT {count} functions; exhaustive p=3,5')
    for p in (3, 5, 7):
        fixed = [s for s in permutations(range(p))
                 if all(s[(i + 1) % p] == (s[i] + 1) % p for i in range(p))]
        translations = {tuple((i + c) % p for i in range(p)) for c in range(p)}
        assert set(fixed) == translations
    print('FIXED_PERMUTATIONS centralizer consists of translations for p=3,5,7')


def main():
    fidelity()
    check_targets()
    check_reduction()
    check_classical_ingredient()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
