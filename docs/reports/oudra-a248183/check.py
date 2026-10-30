#!/usr/bin/env python3
"""Exact finite checks for Oudra's OEIS A248183 conjecture.

Source: https://oeis.org/A248183/internal, revision 19.
Uses only integer arithmetic, Fraction, and integer square roots.
The finite tests are evidence for preregistration, not a universal proof.
"""

from fractions import Fraction
from math import isqrt


SOURCE_TERMS = (
    1, 1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 10, 11, 12,
    12, 13, 14, 15, 15, 16, 17, 17, 18, 19, 20, 20, 21, 22, 22,
    23, 24, 24, 25, 26, 27, 27, 28, 29, 29, 30, 31, 32, 32, 33,
    34, 34, 35, 36, 36, 37, 38, 39, 39, 40, 41, 41, 42, 43, 44,
    44, 45, 46, 46, 47,
)
N_MAX = 10_000


def floor_root_minus_half(q: Fraction) -> int:
    """Compute floor(sqrt(q) - 1/2) using an exact midpoint comparison."""
    assert q >= Fraction(1, 4)
    s = isqrt(q.numerator // q.denominator)
    k = s if Fraction((2 * s + 1) ** 2, 4) <= q else s - 1
    assert Fraction((2 * k + 1) ** 2, 4) <= q
    assert q < Fraction((2 * k + 3) ** 2, 4)
    return k


def main() -> None:
    # This is the literal finite sum, not the proposed telescoped expression.
    remainders = [Fraction(1, 4)]
    for h in range(1, N_MAX + 1):
        remainders.append(remainders[-1] - Fraction(1, h * (h + 1) * (h + 2)))

    # Directly search the source's least positive k; monotonicity permits
    # retaining the previous search position as n increases.
    values = {}
    k = 1
    strict_equalities = []
    for n in range(1, N_MAX + 1):
        threshold = Fraction(1, n * n)
        while remainders[k] >= threshold:
            if remainders[k] == threshold:
                strict_equalities.append((n, k))
            k += 1
        values[n] = k
        assert remainders[k] < threshold
        assert k == 1 or remainders[k - 1] >= threshold

    assert tuple(values[n] for n in range(1, len(SOURCE_TERMS) + 1)) == SOURCE_TERMS
    assert remainders[2] == Fraction(1, 24)
    assert remainders[1] == Fraction(1, 12)
    assert remainders[2] < Fraction(1, 16) < remainders[1]
    assert values[4] == 2
    print(f'FIDELITY: all {len(SOURCE_TERMS)} published DATA terms reproduced by direct finite sums')
    print('FIDELITY: source example a(4)=2; exact inequalities 1/24 < 1/16 < 1/12')

    for h, remainder in enumerate(remainders):
        assert remainder == Fraction(1, 2 * (h + 1) * (h + 2))
    print(f'TAIL: literal finite-sum identity checked for 0 <= k <= {N_MAX}')

    constants = (Fraction(1, 4), Fraction(1, 2), Fraction(1),
                 Fraction(3, 2), Fraction(1749, 1000))
    for n in range(2, N_MAX + 1):
        square_half = Fraction(n * n, 2)
        for c in constants:
            assert values[n] == floor_root_minus_half(square_half + c), (n, c)
    print(f'IDENTITY: Oudra formula checked for every 2 <= n <= {N_MAX}; no discrepancies')
    print('CONSTANTS: same finite range checked at c=1/4, 1/2, 1, 3/2, 1749/1000')

    assert values[12] == 8
    assert remainders[7] == Fraction(1, 144)
    assert floor_root_minus_half(Fraction(12 * 12, 2) + Fraction(24, 100)) == 7
    assert values[3] == 1
    assert floor_root_minus_half(Fraction(3 * 3, 2) + Fraction(7, 4)) == 2
    print('SHARPNESS: n=12 gives 8 at c=1/4 and 7 at c=24/100; n=3 gives 2 at c=7/4 instead of 1')
    print('STRICT_BOUNDARY: equality cases (n, rejected k) = ' + repr(strict_equalities))

    residues = {(2 * j * (j + 1) - n * n) % 4
                for j in range(257) for n in range(257)}
    assert residues == {0, 3}
    print('RESIDUES: 0 <= j,n <= 256 gives exactly {0, 3} modulo 4')
    assert values[1] == 1 and floor_root_minus_half(Fraction(3, 2)) == 0
    print('DOMAIN: n=1 is excluded by the conjecture and its proposed formula gives 0 instead of 1')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
