#!/usr/bin/env python3
"""Exact prospective check of OEIS A183909 (Hardin; Barker).

Source data: https://oeis.org/A183909/internal, revision 8, and its b-file.
OEIS data is attributed to R. H. Hardin and Colin Barker, CC BY-SA 4.0.
This program is experimental evidence, not a Lean proof.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb


MODULUS = 7
CAP = 3
Q = (0, 1, -6, 36, -32, -20, 3, 40, -16, -9, 4)

# R. H. Hardin's displayed OEIS terms, offset 1, revision 8.
DISPLAYED = (
    1, 1, 22, 136, 470, 1193, 2525, 4752, 8238, 13438, 20912,
    31340, 45538, 64475, 89291, 121316, 162090, 213384, 277222,
    355904, 452030, 568525, 708665, 876104, 1074902, 1309554,
    1585020, 1906756, 2280746, 2713535, 3212263, 3784700,
    4439282, 5185148,
)

# Complete table at https://oeis.org/A183909/b183909.txt, n = 1,...,76.
BFILE = (
    1, 1, 22, 136, 470, 1193, 2525, 4752, 8238, 13438,
    20912, 31340, 45538, 64475, 89291, 121316, 162090, 213384, 277222, 355904,
    452030, 568525, 708665, 876104, 1074902, 1309554, 1585020, 1906756, 2280746, 2713535,
    3212263, 3784700, 4439282, 5185148, 6032178, 6991032, 8073190, 9290993, 10657685, 12187456,
    13895486, 15797990, 17912264, 20256732, 22850994, 25715875, 28873475, 32347220, 36161914, 40343792,
    44920574, 49921520, 55377486, 61320981, 67786225, 74809208, 82427750, 90681562, 99612308, 109263668,
    119681402, 130913415, 143009823, 156023020, 170007746, 185021156, 201122890, 218375144, 236842742, 256593209,
    277696845, 300226800, 324259150, 349872974, 377150432, 406176844,
)


def witness_requirements(modulus):
    """Enumerate every unordered two-summand multiset once, plus the target."""
    result = [[] for _ in range(modulus)]
    for left in range(modulus):
        for right in range(left, modulus):
            target = (left + right) % modulus
            needed = [0] * modulus
            needed[target] += 1
            needed[left] += 1
            needed[right] += 1
            result[target].append(tuple(needed))
    return tuple(tuple(row) for row in result)


REQUIREMENTS = witness_requirements(MODULUS)


def valid_multiplicities(counts):
    return all(
        counts[target] == 0
        or any(all(need <= have for need, have in zip(witness, counts))
               for witness in REQUIREMENTS[target])
        for target in range(MODULUS)
    )


def valid_occurrences(values):
    """Independent subset-sum DP, after deleting the represented occurrence.

    No cap or multiplicity-witness predicate is used here. Descending updates
    ensure the two summands are distinct occurrences. Equal target values are
    interchangeable, so it suffices to delete one occurrence of each value.
    """
    seen = set()
    for target_pos, target_value in enumerate(values):
        if target_value in seen:
            continue
        seen.add(target_value)
        reachable = [{0}, set(), set()]
        for pos, value in enumerate(values):
            if pos == target_pos:
                continue
            for taken in (1, 0):
                reachable[taken + 1].update(
                    (residue + value) % MODULUS
                    for residue in tuple(reachable[taken])
                )
        if target_value not in reachable[2]:
            return False
    return True


def formula_coefficient(n):
    """Coefficient of x^n in Q(x)/(1-x)^7, using exact integers."""
    return sum(q * comb(n - degree + 6, 6)
               for degree, q in enumerate(Q) if degree <= n)


def fiber_coefficient(total, histogram):
    result = 0
    for (size, saturated), count in histogram.items():
        if total < size:
            continue
        residual = total - size
        if saturated == 0:
            result += count * int(residual == 0)
        else:
            result += count * comb(residual + saturated - 1, saturated - 1)
    return result


def empirical_polynomial(n):
    return (Fraction(n**6, 720) + Fraction(11 * n**5, 240)
            + Fraction(89 * n**4, 144) + Fraction(209 * n**3, 48)
            - Fraction(13003 * n**2, 360) + Fraction(531 * n, 10) + 12)


def main():
    assert len(DISPLAYED) == 34
    assert BFILE[:34] == DISPLAYED
    assert sum(map(len, REQUIREMENTS)) == comb(MODULUS + 1, 2)
    assert not valid_occurrences((0,))
    assert not valid_occurrences((0, 0))
    assert valid_occurrences((0, 0, 0))

    histogram = Counter()
    by_saturated = Counter()
    numerator = [0] * (CAP * MODULUS + 1)
    states_checked = 0
    valid_small_states = []
    for counts in product(range(CAP + 1), repeat=MODULUS):
        actual_values = tuple(value for value, count in enumerate(counts)
                              for _ in range(count))
        good = valid_multiplicities(counts)
        assert good == valid_occurrences(actual_values), counts
        states_checked += 1
        if not good:
            continue
        size = sum(counts)
        saturated = counts.count(CAP)
        histogram[size, saturated] += 1
        by_saturated[saturated] += 1
        if size <= 2:
            valid_small_states.append(counts)
        for degree in range(MODULUS - saturated + 1):
            numerator[size + degree] += (
                (-1)**degree * comb(MODULUS - saturated, degree))

    target_numerator = [0] * len(numerator)
    for degree in range(MODULUS + 1):
        target_numerator[degree] += (-1)**degree * comb(MODULUS, degree)
    for degree, coefficient in enumerate(Q):
        target_numerator[degree + 2] += coefficient
    assert states_checked == 4**7 == 16384
    assert sum(by_saturated.values()) == 13234
    assert [by_saturated[r] for r in range(8)] == [1239, 3691, 4419, 2729, 945, 189, 21, 1]
    assert valid_small_states == [(0,) * MODULUS]
    assert numerator == target_numerator
    print("CAPPED_STATES", states_checked, "VALID", sum(by_saturated.values()))
    print("VALID_BY_SATURATED", [by_saturated[r] for r in range(8)])
    print("OCCURRENCE_DP_AGREES_ON_ALL_CAPPED_STATES")
    print("FINITE_NUMERATOR_COEFFICIENTS", numerator)
    print("FINITE_IDENTITY_OK: sum_d x^s(d)*(1-x)^(7-r(d)) = (1-x)^7 + x^2*Q(x)")

    direct = []
    arrangements_checked = 0
    n2_solutions = []
    for n in range(1, 11):
        count = 0
        for values in combinations_with_replacement(range(MODULUS), n + 2):
            arrangements_checked += 1
            if valid_occurrences(values):
                count += 1
                if n == 2:
                    n2_solutions.append(values)
        assert count == formula_coefficient(n) == fiber_coefficient(n + 2, histogram), n
        direct.append(count)
    assert n2_solutions == [(0, 0, 0, 0)]
    assert tuple(direct) == DISPLAYED[:10]
    print("SOURCE_EXAMPLE_N2", n2_solutions)
    print("DIRECT_N1_TO_N10", direct)
    print("UNCAPPED_ARRANGEMENTS_CHECKED", arrangements_checked)

    coefficients = [formula_coefficient(n) for n in range(len(BFILE) + 1)]
    assert coefficients[0] == 0
    assert tuple(coefficients[1:35]) == DISPLAYED
    assert tuple(coefficients[1:]) == BFILE
    assert all(fiber_coefficient(n + 2, histogram) == coefficients[n]
               for n in range(1, len(BFILE) + 1))
    assert all(empirical_polynomial(n) == coefficients[n]
               for n in range(4, len(BFILE) + 1))
    for n in range(11, len(BFILE) + 1):
        assert coefficients[n] == sum(
            (-1)**(j + 1) * comb(7, j) * coefficients[n - j]
            for j in range(1, 8))
    onset_residual = coefficients[10] - sum(
        (-1)**(j + 1) * comb(7, j) * coefficients[10 - j]
        for j in range(1, 8))
    assert onset_residual == 4
    print("OEIS_DISPLAYED_TERMS_OK", len(DISPLAYED))
    print("OEIS_BFILE_TERMS_OK", len(BFILE))
    print("FIBER_COEFFICIENTS_AGREE_N1_TO", len(BFILE))
    print("EMPIRICAL_POLYNOMIAL_AGREES_N4_TO", len(BFILE))
    print("RECURRENCE_OK_N11_TO", len(BFILE), "RESIDUAL_AT_N10", onset_residual)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
