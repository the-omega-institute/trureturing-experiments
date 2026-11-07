#!/usr/bin/env python3
"""Exact finite checks for the semantic positive union in OEIS A272266.

Source: Colin Barker, https://oeis.org/A272266, revision 9.
This independent program uses no third-party implementation or recurrence to
construct the sequence. Finite checks are evidence, not an all-index proof.
"""
from fractions import Fraction
from hashlib import sha256
import sys

assert sys.version_info >= (3, 12), sys.version
LIMIT = 30000
SOURCE_DATA = [1, 4, 9, 10, 16, 25, 27, 36, 49, 52, 64, 81, 85, 100, 121, 126, 144, 169, 175, 196, 225, 232, 256, 289, 297, 324, 361, 370, 400, 441, 451, 484, 529, 540, 576, 625, 637, 676, 729, 742, 784, 841, 855, 900, 961, 976, 1024, 1089, 1105, 1156, 1225, 1242, 1296]
BFILE_1000_SHA256 = '02b51581077f19062f21833c0a01faf2eeb5e325a2bde90325a5778009186e67'


def polygonal(m, n):
    """Barker's printed pg formula, with exact rational division."""
    value = Fraction(n * n * (m - 2) - n * (m - 4), 2)
    assert value.denominator == 1
    return value.numerator


def mul(p, q):
    result = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i + j] += a * b
    return result


def power(p, n):
    result = [1]
    for _ in range(n):
        result = mul(result, p)
    return result


def residue_value(n):
    k, r = divmod(n, 3)
    if r == 0:
        assert k >= 1
        return (2 * k + 1) ** 2
    if r == 1:
        return (k + 1) * (4 * k + 1)
    return 4 * (k + 1) ** 2


def main():
    squares = {polygonal(4, n) for n in range(1, LIMIT + 1)}
    decagonals = {polygonal(10, n) for n in range(1, LIMIT + 1)}
    assert all(polygonal(4, n) == n * n for n in range(1, LIMIT + 1))
    assert all(polygonal(10, n) == 4 * n * n - 3 * n for n in range(1, LIMIT + 1))
    # Both omitted tails exceed LIMIT**2. The included squares alone provide
    # LIMIT terms <= LIMIT**2, so sorting these finite sets gives the true prefix.
    assert (LIMIT + 1) ** 2 > LIMIT ** 2
    assert polygonal(10, LIMIT + 1) > LIMIT ** 2
    values = sorted(squares | decagonals)[:LIMIT]
    a = [0] + values
    assert values[:len(SOURCE_DATA)] == SOURCE_DATA
    print(f'FIDELITY: all {len(SOURCE_DATA)} displayed OEIS terms reproduced')
    bfile = ''.join(f'{i} {a[i]}\n' for i in range(1, 1001))
    assert sha256(bfile.encode()).hexdigest() == BFILE_1000_SHA256
    print('FIDELITY: all 1000 b-file rows match their canonical SHA-256')
    assert squares & decagonals == {1}
    assert values.count(1) == 1
    print('OVERLAP: finite square/decagonal sets intersect exactly at 1; union counts it once')
    for k in range(2, LIMIT + 1):
        d = polygonal(10, k)
        assert d - (2 * k - 1) ** 2 == k - 1 > 0
        assert (2 * k) ** 2 - d == 3 * k > 0
    print(f'INTERLEAVING: both strict inequalities checked at k=2..{LIMIT}')
    for n in range(1, LIMIT + 1):
        assert a[n] == residue_value(n), (n, a[n], residue_value(n))
    print(f'ENUMERATION: all three residue formulas match n=1..{LIMIT}')
    for n in range(7, LIMIT + 1):
        assert a[n] - 2 * a[n - 3] + a[n - 6] == 8
    print(f'DIFFERENCE: step-three second difference is 8 at n=7..{LIMIT}')
    for n in range(8, LIMIT + 1):
        assert a[n] == a[n - 1] + 2 * a[n - 3] - 2 * a[n - 4] - a[n - 6] + a[n - 7]
    print(f'RECURRENCE: literal conjecture checked at n=8..{LIMIT}')
    q = mul(power([1, -1], 3), power([1, 1, 1], 2))
    assert q == [1, -1, 0, -2, 2, 0, 1, -1]
    assert q == mul([1, -1], power([1, 0, 0, -1], 2))
    numerator = [0, 1, 3, 5, -1, 0, -1, 1]
    assert numerator == mul([0, 1], [1, 3, 5, -1, 0, -1, 1])
    for n in range(LIMIT + 1):
        coefficient = sum(q[j] * a[n - j] for j in range(min(n, len(q) - 1) + 1))
        assert coefficient == (numerator[n] if n < len(numerator) else 0)
    print('POLYNOMIALS: denominator and numerator factorizations agree exactly in Z[x]')
    print(f'GENERATING FUNCTION: coefficients of QF equal the printed numerator at degrees 0..{LIMIT}')
    assert a[1:9] == [1, 4, 9, 10, 16, 25, 27, 36]
    print('INITIAL: a(1..8) = ' + ','.join(map(str, a[1:9])))
    print('ARITHMETIC: integers and fractions.Fraction only; no floating point')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
