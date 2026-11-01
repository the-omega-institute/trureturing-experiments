#!/usr/bin/env python3
"""Exact finite checks of Hanna's OEIS A326602 square-parity conjecture.

Source: https://oeis.org/A326602, revision 19, FORMULA and displayed data.
No network, floating point, symbolic assumptions, or Lean proof is used.
The finite checks do not establish the unbounded conjecture.
"""

from math import comb, isqrt


DEGREE = 800
SOURCE_TERMS = (
    1, 3, 0, 2, 9, -12, 8, 12, -10, 15, -36, 46, -18, 20, 12, -82,
    101, -162, 302, -168, -100, 92, 32, 40, -244, 351, -452, 1052,
    -1528, 1144, -394, -494, 948, -438, 370, -1474, 2805, -2860,
    2560, -4762, 6554, -4104, 926, 480, -1820, 3074, -2546, 1072,
    -2518, 9745, -17810, 21300, -28982, 37560, -26380, 6162, 686,
    2, 2364, -12342, 30356, -39584, 19448, 7562, 9491, -63824,
    99128, -116668,
)


def reciprocal_one_plus_power(power, last):
    """Coefficients of (1+t)^(-power) by exact polynomial inversion.

    This recurrence is independent of the negative-binomial closed formula
    used in the reindexed computation.
    """
    denominator = [comb(power, h) for h in range(min(power, last) + 1)]
    inverse = [1]
    for k in range(1, last + 1):
        inverse.append(-sum(denominator[h] * inverse[k - h]
                            for h in range(1, min(power, k) + 1)))
    return inverse


def source_expansion(degree):
    """Literal FORMULA (1); n=0 contributes 1/(1+x)."""
    result = reciprocal_one_plus_power(1, degree)
    for n in range(1, degree + 1):
        p = 2 * n + 1
        inverse = reciprocal_one_plus_power(n + 1, (degree - n) // p)
        for j in range(min(n, (degree - n) // (2 * n - 1)) + 1):
            base = n + (2 * n - 1) * j
            numerator = comb(n, j)
            for k in range((degree - base) // p + 1):
                result[base + p * k] += numerator * inverse[k]
    return result


def symmetric_expansion(degree):
    """Independently enumerate all triples using exact multinomial weights."""
    result = [0] * (degree + 1)
    triples = 0
    for j in range(isqrt(degree // 2) + 1):
        for l in range(degree + 1):
            base = l + 2 * (l + j) * j
            if base > degree:
                break
            for k in range((degree - base) // (1 + 2 * (l + j)) + 1):
                d = l + k + 2 * (l + j) * (j + k)
                weight = comb(l + j + k, k) * comb(l + j, j)
                swapped_d = k + l + 2 * (k + j) * (j + l)
                swapped_weight = comb(k + j + l, l) * comb(k + j, j)
                assert d == swapped_d
                assert weight == swapped_weight
                assert d >= l + k + 2 * j * j
                assert d == (l + j) + (2 * (l + j) - 1) * j + (2 * (l + j) + 1) * k
                if l == k and l > 0:
                    assert weight == comb(2 * l + j, j) * comb(2 * l, l)
                    assert comb(2 * l, l) == 2 * comb(2 * l - 1, l - 1)
                    assert weight % 2 == 0
                if d % 2:
                    assert l != k and (l + k) % 2 == 1
                    assert (-1) ** k * weight + (-1) ** l * swapped_weight == 0
                result[d] += (-1) ** k * weight
                triples += 1
    return result, triples


def main():
    literal = source_expansion(DEGREE)
    symmetric, triples = symmetric_expansion(DEGREE)
    assert literal == symmetric
    coefficients = literal[::2]
    assert tuple(coefficients[:len(SOURCE_TERMS)]) == SOURCE_TERMS
    assert all(literal[d] == 0 for d in range(1, DEGREE + 1, 2))
    assert all((value % 2 == 1) == (isqrt(m) ** 2 == m)
               for m, value in enumerate(coefficients))
    assert all((value % 2 == 1) == (d % 2 == 0 and isqrt(d // 2) ** 2 == d // 2)
               for d, value in enumerate(literal))
    print(f"FIDELITY: all {len(SOURCE_TERMS)} displayed OEIS terms match")
    print(f"SOURCE_VS_SYMMETRIC: all degrees 0..{DEGREE} match; {triples} triples")
    print(f"ODD_DEGREES: all {DEGREE // 2} coefficients vanish")
    print(f"SQUARE_PARITY: a(0)..a({DEGREE // 2}) all match")
    print("ODD_INDICES:", [m for m, value in enumerate(coefficients) if value % 2])
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
