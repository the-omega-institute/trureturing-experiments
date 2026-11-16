#!/usr/bin/env python3
"""Exact finite checks of OEIS A316553; independently written, no CAS needed.

Source: https://oeis.org/A316553, revision 23; see README.md for attribution.
No floating point, CRT multiplication, or conjectured local count is used in
constructing the matrix sets. Full four-entry searches cross-check pruning.
"""
from itertools import product
from math import gcd
from pathlib import Path


def square(matrix: tuple[int, int, int, int], n: int) -> tuple[int, ...]:
    a, b, c, d = matrix
    return ((a*a+b*c) % n, (a*b+b*d) % n,
            (c*a+d*c) % n, (c*b+d*d) % n)


def is_square_one_sl2(matrix: tuple[int, int, int, int], n: int) -> bool:
    a, b, c, d = matrix
    return (a*d-b*c) % n == 1 % n and square(matrix, n) == (1 % n, 0, 0, 1 % n)


def all_four_entries(n: int) -> set[tuple[int, int, int, int]]:
    return {m for m in product(range(n), repeat=4) if is_square_one_sl2(m, n)}


def adjugate_pruned(n: int) -> set[tuple[int, int, int, int]]:
    # A in SL2 has inverse adj(A). A^2=I thus forces A=adj(A):
    # d=a, b=-b and c=-c. Enumerate this entire necessary domain,
    # then independently test the original determinant and matrix square.
    ann2 = [x for x in range(n) if (2*x) % n == 0]
    return {(a, b, c, a) for a in range(n) for b in ann2 for c in ann2
            if is_square_one_sl2((a, b, c, a), n)}


def conjectured_a(n: int) -> int:
    remainder, p, distinct, v2 = n, 2, 0, 0
    while p*p <= remainder:
        if remainder % p == 0:
            distinct += 1
            exponent = 0
            while remainder % p == 0:
                remainder //= p
                exponent += 1
            if p == 2:
                v2 = exponent
        p += 1
    if remainder > 1:
        distinct += 1
        if remainder == 2:
            v2 = 1
    return 2**(distinct + min(3, v2)) - 1


def main() -> None:
    # Fidelity: the source's three explicitly displayed matrices at n=2.
    displayed = {(0, 1, 1, 0), (1, 0, 1, 1), (1, 1, 0, 1)}
    assert all_four_entries(2) - {(1, 0, 0, 1)} == displayed
    assert all_four_entries(1) == {(0, 0, 0, 0)}
    assert (1 % 2, 0, 0, 1 % 2) == (-1 % 2, 0, 0, -1 % 2)
    print('FIDELITY n=1: a(1)=0; n=2: exactly the three displayed matrices; -I=I')

    # Check the adjugate characterization for every raw matrix at small n.
    for n in range(1, 21):
        for a, b, c, d in product(range(n), repeat=4):
            lhs = is_square_one_sl2((a, b, c, d), n)
            rhs = (d == a and 2*b % n == 0 and 2*c % n == 0
                   and (a*a-b*c) % n == 1 % n)
            assert lhs == rhs, (n, a, b, c, d, lhs, rhs)
    print('CHARACTERIZATION all n=1..20, every one of n^4 matrices: OK')

    for n in range(1, 33):
        assert all_four_entries(n) == adjugate_pruned(n), n
    print('UNRESTRICTED_CROSS_CHECK all n=1..32: identical matrix sets')

    counts = {}
    for n in range(1, 2001):
        matrices = adjugate_pruned(n)
        identity = (1 % n, 0, 0, 1 % n)
        assert identity in matrices
        counts[n] = len(matrices - {identity})
        assert counts[n] == conjectured_a(n), (n, counts[n], conjectured_a(n))
    print('CLAIM exhaustive adjugate-pruned enumeration, all n=1..2000: OK')

    source_rows = []
    for line in Path(__file__).with_name('oeis-b316553.txt').read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            n, value = map(int, line.split())
            source_rows.append((n, value))
            assert counts[n] == value, (n, counts[n], value)
    assert [n for n, _ in source_rows] == list(range(1, 1024))
    print(f'OEIS_TABLE all {len(source_rows)} rows n=1..1023: exact match')

    values = []
    for e in range(1, 9):
        n = 2**e
        matrices = adjugate_pruned(n)
        assert len(matrices) == 2**(1 + min(3, e))
        if e >= 2:
            assert all((b*c) % n == 0 for _, b, c, _ in matrices)
            roots = [a for a in range(n) if a*a % n == 1]
            assert len(roots) == (2 if e == 2 else 4)
            if e >= 3:
                lifts = {a for a in range(n)
                         if a % (2**(e-1)) in {1, 2**(e-1)-1}}
                assert set(roots) == lifts
        values.append(len(matrices))
    print('POWERS_OF_TWO N(2^e), e=1..8:', ','.join(map(str, values)))

    pairs = 0
    for m in range(1, 51):
        for n in range(1, 51):
            if m*n <= 2000 and gcd(m, n) == 1:
                assert counts[m*n]+1 == (counts[m]+1)*(counts[n]+1)
                pairs += 1
    print(f'ISRAEL_MULTIPLICATIVITY {pairs} coprime ordered pairs m,n<=50, mn<=2000: OK')
    print('COUNTS n,a(n),N(n) for all checked moduli')
    for n in range(1, 2001):
        print(f'{n} {counts[n]} {counts[n]+1}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
