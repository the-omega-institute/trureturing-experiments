#!/usr/bin/env python3
"""Exact prospective checks for Hardin--Barker OEIS A183901.

Source: https://oeis.org/A183901/internal, revision 8.
All arithmetic used as evidence is integer arithmetic. No third-party packages.
"""
from itertools import combinations_with_replacement, product
from math import comb

MODULUS = 7
CAP = 4
# Ascending coefficients of the whole numerator, including its leading x.
Q = (0, 1, -6, 150, -302, -72, 649, -548, 60, 102, -33)
# Hardin's displayed DATA, offset 1, and four further terms from his b-file.
DISPLAYED = (
    1, 1, 136, 664, 1720, 3491, 6263, 10400, 16357, 24694, 36091,
    51364, 71482, 97585, 131003, 173276, 226175, 291724, 372223,
    470272, 588796, 731071, 900751, 1101896, 1339001, 1617026,
    1941427, 2318188, 2753854, 3255565, 3831091, 4488868, 5238035,
)
BFILE_TAIL = (6088472, 7050839, 8136616, 9358144)


def witness_demands():
    """All 84 size-three value multisets, grouped by their sum modulo 7.

    Add one target occurrence to the three witness occurrences. Keeping only
    positive coordinates is an exact sparse representation of each demand.
    """
    demands = [[] for _ in range(MODULUS)]
    for triple in combinations_with_replacement(range(MODULUS), 3):
        target = sum(triple) % MODULUS
        need = [0] * MODULUS
        need[target] += 1
        for value in triple:
            need[value] += 1
        demands[target].append(tuple((j, v) for j, v in enumerate(need) if v))
    assert sum(map(len, demands)) == comb(MODULUS + 2, 3) == 84
    return demands


DEMANDS = witness_demands()


def capped_valid(c):
    return all(
        c[i] == 0 or any(all(c[j] >= v for j, v in need) for need in DEMANDS[i])
        for i in range(MODULUS)
    )


def occurrence_valid(arrangement):
    """Independent literal check using distinct individual positions.

    Remove a target position, then compute reachable residues after selecting
    exactly three of the remaining positions. A bit represents one residue.
    Descending updates prevent selecting one position more than once.
    """
    for target in set(arrangement):
        removed = False
        reach = [1, 0, 0, 0]
        for value in arrangement:
            if value == target and not removed:
                removed = True
                continue
            for k in (3, 2, 1):
                bits = reach[k - 1]
                rotated = ((bits << value) | (bits >> (MODULUS - value))) & 127
                reach[k] |= rotated
        if not reach[3] & (1 << target):
            return False
    return True


def source_coefficient(n):
    return sum(Q[k] * comb(n - k + 6, 6) for k in range(min(n, 10) + 1))


def main():
    valid = []
    by_r = [0] * 8
    polynomial = [0] * 29
    for d in product(range(CAP + 1), repeat=MODULUS):
        if not capped_valid(d):
            continue
        s = sum(d)
        r = d.count(CAP)
        valid.append((s, r))
        by_r[r] += 1
        for k in range(8 - r):
            polynomial[s + k] += (-1) ** k * comb(7 - r, k)
    assert len(valid) == 74311
    assert by_r == [14267, 27206, 21273, 8960, 2240, 336, 28, 1]
    expected = [0] * 29
    for k in range(8):
        expected[k] = (-1) ** k * comb(7, k)
    for k, coefficient in enumerate(Q):
        expected[k + 3] += coefficient
    assert polynomial == expected
    assert [(s, r) for s, r in valid if s <= 3] == [(0, 0)]
    print('capped_states_checked:', 5 ** 7)
    print('valid_capped_states:', len(valid))
    print('valid_counts_by_saturated_coordinates:', by_r)
    print('finite_identity_coefficients_degrees_0_through_28:', polynomial)

    independent = []
    checked = []
    for n in range(1, 8):
        count = 0
        number = 0
        for a in combinations_with_replacement(range(7), n + 3):
            number += 1
            literal = occurrence_valid(a)
            c = tuple(a.count(i) for i in range(7))
            assert literal == capped_valid(c) == capped_valid(tuple(min(v, CAP) for v in c))
            count += literal
        assert number == comb(n + 9, 6)
        assert count == DISPLAYED[n - 1] == source_coefficient(n)
        independent.append(count)
        checked.append(number)
    print('independent_occurrence_checks_n_1_through_7:', independent)
    print('sorted_tuples_checked_by_n:', checked)
    print('source_fidelity_a_3:', independent[2])

    coefficients = [source_coefficient(n) for n in range(1, 38)]
    assert tuple(coefficients[:33]) == DISPLAYED
    assert tuple(coefficients[33:]) == BFILE_TAIL
    # Compute coefficients of the cap-fiber series independently of Q.
    for n in range(1, 38):
        size = n + 3
        from_fibers = sum(
            (int(size == s) if r == 0 else comb(size - s + r - 1, r - 1))
            for s, r in valid if s <= size
        )
        assert from_fibers == coefficients[n - 1]
    print('source_bfile_and_cap_fibers_n_1_through_37:', coefficients)
    recurrence_weights = (7, -21, 35, -35, 21, -7, 1)
    for n in range(11, 301):
        assert source_coefficient(n) == sum(
            w * source_coefficient(n - j)
            for j, w in enumerate(recurrence_weights, 1)
        )
    # n=10 is excluded in the source, and is indeed not a homogeneous case.
    residual_10 = source_coefficient(10) - sum(
        w * source_coefficient(10 - j)
        for j, w in enumerate(recurrence_weights, 1)
    )
    assert residual_10 == -33
    print('recurrence_checked_n_11_through_300: OK')
    print('recurrence_residual_at_n_10:', residual_10)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
