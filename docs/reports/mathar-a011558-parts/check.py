#!/usr/bin/env python3
"""Exact finite checks of Mathar's OEIS A011558 partition conjecture.

Source: https://oeis.org/A011558 (revision 76), https://oeis.org/A084215
(revision 51). Coefficient lists below are source DATA, starting at offset 0.
Uses only Python integers. These checks do not prove the unbounded conjecture.
"""
from collections import Counter
from itertools import combinations

LIMIT = 10000
SOURCE_A011558 = [0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0]
SOURCE_A084215 = [1, 2, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840, 327680, 655360, 1310720, 2621440, 5242880, 10485760, 20971520, 41943040, 83886080, 167772160, 335544320, 671088640, 1342177280, 2684354560, 5368709120, 10737418240]


def source_coefficients(last):
    # (1+x+x^2+x^3+x^4)*A = x+x^3 over F_2.
    a = []
    for n in range(last + 1):
        numerator = int(n in (1, 3))
        a.append((numerator - sum(a[n-j] for j in range(1, min(n, 4)+1))) % 2)
    return a


def part_coefficients(last):
    # (1-2*x)*B = 1+x^2 over Z, without the conjectured closed form.
    b = []
    for n in range(last + 1):
        b.append(int(n in (0, 2)) + (2*b[n-1] if n else 0))
    return b


def main():
    a = source_coefficients(LIMIT + 1)
    b = part_coefficients(max(40, len(SOURCE_A084215)-1))
    assert a[:len(SOURCE_A011558)] == SOURCE_A011558
    assert b[:len(SOURCE_A084215)] == SOURCE_A084215
    assert (len(SOURCE_A011558), len(SOURCE_A084215)) == (106, 34)
    assert b[:3] == [1, 2, 5]
    assert all(b[n] == 2*b[n-1] for n in range(3, len(b)))
    assert all(b[k] == 5*2**(k-2) for k in range(2, len(b)))
    print('Fidelity: A011558 DATA 106 terms; A084215 DATA 34 terms: OK')
    print('Fidelity: A084215 stated doubling recurrence and initial values: OK')

    # Parts come from source B's coefficient recurrence, not from n mod 5.
    parts = sorted(set(t for t in b if t <= LIMIT))
    assert len(parts) == len([t for t in b if t <= LIMIT])
    assert parts == [1, 2] + [5*2**j for j in range(40) if 5*2**j <= LIMIT]
    dp = [0]*(LIMIT+1)
    dp[0] = 1
    for t in parts:
        # Descending order counts each value at most once.
        for n in range(LIMIT, t-1, -1):
            dp[n] += dp[n-t]

    # Independently enumerate every finite subset of the complete bounded part set.
    # Positive parts greater than LIMIT cannot occur in a partition of n <= LIMIT.
    enumeration = Counter()
    subsets = 0
    for size in range(len(parts)+1):
        for subset in combinations(parts, size):
            subsets += 1
            enumeration[sum(subset)] += 1
    assert subsets == 2**len(parts)
    assert all(dp[n] == enumeration[n] for n in range(LIMIT+1))
    assert all(dp[n] == a[n+1] for n in range(LIMIT+1))
    assert all(dp[n] == int(n % 5 != 4) for n in range(LIMIT+1))
    assert dp[0] == 1  # empty partition is included
    print('Allowed parts <= 10000:', parts)
    print('All subsets enumerated:', subsets)
    print('Identity: DP = exhaustive subsets = original F_2 coefficient a(n+1), 0 <= n <= 10000: OK')
    print('Residue classification: zero on residue 4; one on residues 0,1,2,3: OK')
    print('Counts on 0..10000:', dict(sorted(Counter(dp).items())))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
