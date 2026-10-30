#!/usr/bin/env python3
"""Exact independent finite checks for Paul D. Hanna's OEIS A292180.

The source FORMULA is expanded directly by source_coefficient. A second
implementation multiplies truncated copies of (1+t)/(1-t), whose positive
coefficients are all 2. No conjectural divisibility or parity is assumed.
"""
from hashlib import sha256
from math import comb, isqrt
from pathlib import Path

LIMIT = 1024
SOURCE_FIRST_40 = (
    4, 0, 16, 16, 24, 0, 32, 96, 116, 0, 48, 192, 56, 0, 608, 704,
    72, 0, 80, 480, 1408, 0, 96, 3712, 2108, 0, 2720, 896, 120, 0,
    128, 9600, 4672, 0, 17088, 12112, 152, 0, 7392, 20800,
)
SOURCE_SQUARE_QUOTIENTS = (
    1, 4, 29, 176, 527, 3028, 14457, 107200, 446745, 2392604,
    13286165, 140564336, 415382567, 2333455268, 17078911507,
    78663453440, 419472490547, 2377516612900, 13482186743565,
    78663154105296, 437169506932981, 2481447593907572,
    14146164790774889, 161511806183206336, 460995825168188653,
    2634869356953946428, 15071070681878977525, 86632929673574593072,
    494051395886263605335, 2955861929786748934348,
    16234283204352299108321,
)


def source_coefficient(d, k):
    """[t^k] ((1+t)^(2d)+(-1)^d(1-t)^(2d))/(1-t^2)^d."""
    numerator_times_inverse = sum(
        comb(d + ell - 1, ell) * comb(2 * d, k - 2 * ell)
        for ell in range(k // 2 + 1)
        if k - 2 * ell <= 2 * d
    )
    return (1 + (-1) ** (d + k)) * numerator_times_inverse


def direct_source_coefficients(limit):
    result = [0] * (limit + 1)
    for d in range(1, limit + 1):
        for k in range(1, limit // d + 1):
            result[d * k] += source_coefficient(d, k)
    return result


def reflected_product_coefficients(limit):
    """Independent coefficient convolution for F_d and its reflected copy."""
    previous = [1] + [0] * limit
    result = [0] * (limit + 1)
    checked = 0
    for d in range(1, limit + 1):
        b = []
        prefix = 0
        for k in range(limit // d + 1):
            b.append(previous[k] + 2 * prefix)
            prefix += previous[k]
        for k in range(1, len(b)):
            assert (b[k] - 2 * d) % 4 == 0, (d, k, b[k])
            checked += 1
            result[d * k] += (1 + (-1) ** (d + k)) * b[k]
        previous = b
    return result, checked


def main():
    source_file = Path(__file__).with_name('b292180.txt')
    raw = source_file.read_bytes()
    table = {}
    for line in raw.decode('utf-8-sig').splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            n, value = map(int, line.split())
            assert n not in table
            table[n] = value
    assert sorted(table) == list(range(1, LIMIT + 1))

    direct = direct_source_coefficients(LIMIT)
    reflected, local_checks = reflected_product_coefficients(LIMIT)
    assert direct == reflected, 'independent coefficient algorithms disagree'
    assert tuple(direct[1:41]) == SOURCE_FIRST_40
    assert all(direct[n] == table[n] for n in table)
    for root, q in enumerate(SOURCE_SQUARE_QUOTIENTS, 1):
        assert direct[root * root] == 4 * q, (root, q)
    print('Source: Paul D. Hanna, OEIS A292180, revision 15; offset 1.')
    print('b-file sha256:', sha256(raw).hexdigest())
    print('Fidelity: all 1024 b-file values, first 40 displayed coefficients,')
    print('and all 31 published square quotients agree exactly.')
    print('Independent FORMULA expansion and reflected-product convolution agree at N=0..1024.')
    print('Coefficient congruence b_d(k) = 2d (mod 4):', local_checks, 'positive (d,k) pairs with d*k<=1024.')

    divisors = [[] for _ in range(LIMIT + 1)]
    for d in range(1, LIMIT + 1):
        for n in range(d, LIMIT + 1, d):
            divisors[n].append(d)
    for n in range(1, LIMIT + 1):
        assert direct[n] % 4 == 0
        if n % 4 == 2:
            assert direct[n] == 0
        if n % 2 == 1:
            sigma = sum(divisors[n])
            assert (direct[n] - 4 * sigma) % 8 == 0
            assert (direct[n] // 4) % 2 == (isqrt(n) ** 2 == n)
    odd_roots = list(range(1, isqrt(LIMIT) + 1, 2))
    for root in odd_roots:
        assert direct[root * root] % 8 == 4
    print('All 1024 coefficients are divisible by 4; all N=2 (mod 4) coefficients vanish.')
    print('All 512 positive odd N<=1024 satisfy a(N)=4*sigma_1(N) (mod 8)')
    print('and a(N)/4 is odd exactly when N is a square.')
    print('Published conjecture: all', len(odd_roots), 'odd squares <=1024 pass.')
    print('Odd-root quotients for roots 1,3,5,7,9:', [direct[r * r] // 4 for r in (1, 3, 5, 7, 9)])
    print('These are exact finite checks, not a proof for all indices.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
