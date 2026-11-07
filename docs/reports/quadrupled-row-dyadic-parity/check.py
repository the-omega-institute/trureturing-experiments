#!/usr/bin/env python3
"""Exact finite checks for OEIS A397594 (y=4)."""
from fractions import Fraction

SOURCE_PREFIX = [
    1,
    3,
    27,
    528,
    16512,
    705195,
    37887834,
    2442170091,
    183244178880,
    15675764998821,
    1505681872695420,
    160476900178362804,
    18799367838421680984,
    2401670597893104602064,
    332382395344108056499278,
    49549320572895140393543289,
    7916911857698298814547857788,
    1349889314907247276753745833113,
    244674159394172002450191865114968,
    46982018676209283983413873901242380,
    9527740460353118279837118091928699196,
    2034963272634828858009326663626862180274,
    456604045240145780477191741291813105800888,
    107385782328038188275949085412202948700339886,
    26416274035545447424931364352419708827660220132,
]


def mul(a, b, limit):
    out = [0] * (min(limit, len(a) + len(b) - 2) + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            d = i + j
            if d > limit:
                break
            out[d] += ai * bj
    return out


def power(a, exponent, limit):
    result = [1]
    base = a[: limit + 1]
    while exponent:
        if exponent & 1:
            result = mul(result, base, limit)
        exponent >>= 1
        if exponent:
            base = mul(base, base, limit)
    return result


def row(a, n):
    """[x^n] A(x)^(n+1)/(1-(n+1)x), exactly over Z."""
    m = n + 1
    p = power(a[: n + 1], m, n)
    return sum(p[j] * m ** (n - j) for j in range(n + 1))


def normalized_source_sum(a, n):
    """Sum_{j=0}^n [x^j] A(x/(n+1))^(n+1), exactly over Q."""
    m = n + 1
    p = power(a[: n + 1], m, n)
    return sum(Fraction(p[j], m**j) for j in range(n + 1))


def construct(count):
    coefficients = [1]
    for n in range(1, count):
        m = n + 1
        # The new coefficient is set to zero while computing the known part.
        known = row(coefficients + [0], n)
        target = (4 * m) ** n
        numerator = target - known
        assert numerator % m == 0, (n, numerator, m)
        coefficients.append(numerator // m)
    return coefficients


def expected_odd(n):
    return any(n == (1 << k) - 1 or n == (1 << k) + 1 for k in range(2, n.bit_length() + 2))


def main():
    limit = 101
    coefficients = construct(limit)

    # Source fidelity: all 25 values in the OEIS b-file prefix agree.
    assert coefficients[: len(SOURCE_PREFIX)] == SOURCE_PREFIX
    print("OEIS A397594 b-file prefix reproduced: degrees 0..24")

    # Both equivalent source normalizations hold at every tested degree.
    for n in range(limit - 1):
        assert row(coefficients, n) == (4 * (n + 1)) ** n
    for n in range(40):
        assert normalized_source_sum(coefficients, n) == Fraction(4**n)
    print("exact source identities checked: formula (1) through n=39; formula (2) through n=100")

    # The proposed dyadic-neighbor support has no finite counterexample in 4..100.
    odd_indices = [n for n in range(4, limit) if coefficients[n] & 1]
    expected_indices = [n for n in range(4, limit) if expected_odd(n)]
    assert odd_indices == expected_indices
    print("odd indices for 4 <= n <= 100:", odd_indices)
    print("first coefficients:", coefficients[:10])
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
