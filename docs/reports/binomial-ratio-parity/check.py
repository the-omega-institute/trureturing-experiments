#!/usr/bin/env python3
"""Exact finite checks for Peter Bala's OEIS A295462 parity conjecture.

The checks supply finite evidence, not an all-n proof. No third-party modules.
"""
from fractions import Fraction
from math import factorial, prod

SOURCE_TERMS = (
    1, 1001, 2731365, 8452632980, 27614258931525,
    92994778322102376, 319323042934810327908,
    1111423634892509277158820, 3907116913932203786698907205,
    13840531451943824412017357524400,
    49326067264186310295180241747943240,
    176657291645522807951014259080005974280,
)
TABLE = (0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
         0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1)
RATIO_LIMIT = 4095
VALUATION_LIMIT = 65535


def ratio(n):
    return Fraction(factorial(15*n) * factorial(2*n),
                    factorial(10*n) * factorial(6*n) * factorial(n))


def spaced(n):
    positions = [i for i in range(n.bit_length()) if (n >> i) & 1]
    return all(t-e >= 4 for e, t in zip(positions, positions[1:]))


def floor_defect(n, q):
    return 15*n//q + 2*n//q - 10*n//q - 6*n//q - n//q


def landau(x):
    return sum(c * ((a*x).numerator // (a*x).denominator)
               for a, c in ((15, 1), (2, 1), (10, -1), (6, -1), (1, -1)))


def vp_factorial(n, p):
    total = 0
    while n:
        n //= p
        total += n
    return total


def v2_ratio(n):
    return (vp_factorial(15*n, 2) + vp_factorial(2*n, 2)
            - vp_factorial(10*n, 2) - vp_factorial(6*n, 2)
            - vp_factorial(n, 2))


def direct_v2(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def main():
    # Fidelity: OEIS %S/%T/%U initial values and Bala's proved divisibilities.
    assert tuple(ratio(n) for n in range(len(SOURCE_TERMS))) == SOURCE_TERMS
    print("Fidelity: all 12 published A295462 terms reproduced exactly (n=0..11).")

    # Every discontinuity of F in [0,1) occurs on this grid. Thus the value
    # at a left endpoint determines the whole half-open cell.
    breaks = {Fraction(j, a) for a in (15, 2, 10, 6, 1) for j in range(a)}
    assert all((30*x).denominator == 1 for x in breaks)
    actual = tuple(landau(Fraction(j, 30)) for j in range(30))
    assert actual == TABLE and set(actual) == {0, 1}
    assert all(landau(Fraction(2*j+1, 60)) == TABLE[j] for j in range(30))
    assert all(landau(Fraction(j, 30)+k) == TABLE[j]
               for j in range(30) for k in (-3, -1, 0, 1, 4))
    print("Landau 30-cell table:", list(actual))

    # Check all endpoint floor bounds using Fractions and strict upper ends.
    intervals = ((1, Fraction(3, 4), Fraction(23, 30), (11, 1, 7, 4)),
                 (2, Fraction(5, 8), Fraction(19, 30), (9, 1, 6, 3)),
                 (3, Fraction(9, 16), Fraction(17, 30), (8, 1, 5, 3)))
    for d, lo, hi, expected in intervals:
        at_lo = tuple((a*lo).numerator // (a*lo).denominator for a in (15, 2, 10, 6))
        below_hi = tuple(((a*hi).numerator-1) // (a*hi).denominator
                         for a in (15, 2, 10, 6))
        assert at_lo == below_hi == expected
        assert sum((1, 1, -1, -1)[i]*expected[i] for i in range(4)) == 1
        print(f"First bad gap {d}: interval [{lo}, {hi}), floors {expected}, F=1.")

    value = 1
    odd_count = 0
    for n in range(RATIO_LIMIT+1):
        if n:
            numerator = prod(range(15*(n-1)+1, 15*n+1))
            numerator *= prod(range(2*(n-1)+1, 2*n+1))
            denominator = prod(range(10*(n-1)+1, 10*n+1))
            denominator *= prod(range(6*(n-1)+1, 6*n+1)) * n
            value, remainder = divmod(value*numerator, denominator)
            assert remainder == 0
        if n <= 256:
            assert ratio(n) == value
        assert value % (6*n+1) == 0 and value % (10*n+1) == 0
        assert bool(value % 2) == spaced(n)
        assert direct_v2(value) == v2_ratio(n)
        odd_count += value % 2
    print(f"Exact integer recurrence: integrality, both source divisibilities, parity,")
    print(f"and direct v2 agree for all n=0..{RATIO_LIMIT}; odd indices={odd_count}.")
    print("Independent factorial/Fraction evaluation agrees for n=0..256.")

    witnesses = [0, 0, 0]
    sparse_count = 0
    largest_v2 = 0
    for n in range(VALUATION_LIMIT+1):
        valuation = v2_ratio(n)
        digit_value = (5*n).bit_count() + (3*n).bit_count() - (15*n).bit_count()
        terms = [floor_defect(n, 1 << k) for k in range(1, (15*n).bit_length()+1)]
        assert valuation == digit_value == sum(terms)
        assert all(v in (0, 1) for v in terms)
        assert (valuation == 0) == spaced(n)
        largest_v2 = max(largest_v2, valuation)
        positions = [i for i in range(n.bit_length()) if (n >> i) & 1]
        bad = [(e, t) for e, t in zip(positions, positions[1:]) if t-e < 4]
        if bad:
            e, t = bad[0]
            d = t-e
            lower = n % (1 << e)
            assert 15*lower < (1 << e)
            x = Fraction(n % (1 << (t+1)), 1 << (t+1))
            _, lo, hi, _ = intervals[d-1]
            assert lo <= x < hi and landau(x) == floor_defect(n, 1 << (t+1)) == 1
            assert valuation >= 1
            witnesses[d-1] += 1
        else:
            sparse_count += 1
            r = n.bit_count()
            assert (3*n).bit_count() == (5*n).bit_count() == 2*r
            assert (15*n).bit_count() == 4*r
    print(f"Valuation/digit-sum/spacing checks: all n=0..{VALUATION_LIMIT};")
    print(f"sparse={sparse_count}, first-bad gap counts={witnesses}, max v2={largest_v2}.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
