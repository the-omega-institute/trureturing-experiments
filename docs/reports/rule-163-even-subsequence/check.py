#!/usr/bin/env python3
"""Exact finite checks of the OEIS Rule 163 even-subsequence conjecture.

No dependencies, no network access, and no floating-point arithmetic.
Finite corroboration only: this program does not establish an all-time theorem.
"""

from fractions import Fraction


LIMIT = 300
OEIS_A266753 = (
    1, 4, 18, 74, 298, 1194, 4778, 19114, 76458, 305834, 1223338,
    4893354, 19573418, 78293674, 313174698, 1252698794, 5010795178,
    20043180714, 80172722858, 320690891434, 1282763565738,
    5131054262954, 20524217051818, 82096868207274, 328387472829098,
)
OEIS_A392003 = (
    0, 4, 7, 18, 35, 74, 147, 298, 595, 1194, 2387, 4778,
    9555, 19114, 38227, 76458, 152915, 305834, 611667, 1223338,
    2446675, 4893354, 9786707, 19573418, 39146835, 78293674,
    156587347, 313174698, 626349395, 1252698794, 2505397587,
    5010795178, 10021590355, 20043180714,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def table_update(left, centre, right):
    return (163 >> (4 * left + 2 * centre + right)) & 1


def boolean_update(left, centre, right):
    # Independent expression for the truth table, not a row formula.
    return int((left and right) or (not left and not centre))


def B(m):
    return Fraction(7 * 2**m - 9 + 5 * (-1)**m, 6)


def run():
    table = tuple(table_update((k >> 2) & 1, (k >> 1) & 1, k & 1)
                  for k in range(7, -1, -1))
    require(table == (1, 0, 1, 0, 0, 0, 1, 1), 'Wolfram rule numbering')
    for k in range(8):
        bits = ((k >> 2) & 1, (k >> 1) & 1, k & 1)
        require(table_update(*bits) == boolean_update(*bits), 'truth-table cross-check')

    # Simulator 1: exact finite exceptions from the spatially uniform background.
    # Its background evolves 0 -> 1 -> 1; it never imposes a zero exterior.
    background = 0
    exceptions = {0: 1}

    # Simulator 2: a padded dense row. Shrink by one on both ends after each step.
    # Thus every retained cell has its complete causal cone in the original row.
    radius = 2 * LIMIT + 4
    dense = [int(i == 0) for i in range(-radius, radius + 1)]
    lo = -radius
    values = []
    cell_checks = 0
    for n in range(LIMIT + 1):
        def cell(i):
            return exceptions.get(i, background)

        for i in range(-n - 3, n + 4):
            require(cell(i) == dense[i - lo], f'simulators disagree: n={n}, i={i}')
            cell_checks += 1
        word = [cell(i) for i in range(-n, n + 1)]
        value = 0
        for bit in word:
            value = 2 * value + bit
        require(value == sum(cell(n - j) * 2**j for j in range(2 * n + 1)),
                f'encoding orientation n={n}')
        values.append(value)
        if n >= 1:
            require(word == [1, 0, 0] + [1, 0] * (n - 1), f'row word n={n}')
            require(background == 1, f'nonquiescent exterior n={n}')
            require(all(-n <= i <= n for i in exceptions), f'exterior support n={n}')
            require(B(2 * n) == value, f'full subsequence identity n={n}')
            require(Fraction(7 * 4**n - 4, 6) == value, f'closed form n={n}')
            require(value == 4**n + 2 * sum(4**r for r in range(n - 1)),
                    f'geometric sum n={n}')
        if n >= 2:
            require(value == 4 * values[n - 1] + 2, f'first recurrence n={n}')
        if n >= 3:
            require(value == 5 * values[n - 1] - 4 * values[n - 2],
                    f'second recurrence n={n}')
        if n < LIMIT:
            candidates = {i + d for i in exceptions for d in (-1, 0, 1)}
            new_background = table_update(background, background, background)
            new_exceptions = {}
            for i in candidates:
                bit = table_update(cell(i - 1), cell(i), cell(i + 1))
                if bit != new_background:
                    new_exceptions[i] = bit
            background, exceptions = new_background, new_exceptions
            dense = [boolean_update(dense[j - 1], dense[j], dense[j + 1])
                     for j in range(1, len(dense) - 1)]
            lo += 1

    require(tuple(values[:len(OEIS_A266753)]) == OEIS_A266753, 'source A266753 data')
    require(tuple(B(m) for m in range(1, 35)) == OEIS_A392003, 'source A392003 data')
    # (1-X)(1-4X) * sum A(n)X^n = 1-X+2X^2, coefficient by coefficient.
    for n, value in enumerate(values):
        coeff = value - (5 * values[n - 1] if n >= 1 else 0)
        coeff += 4 * values[n - 2] if n >= 2 else 0
        require(coeff == {0: 1, 1: -1, 2: 2}.get(n, 0), f'GF coefficient n={n}')

    print('Rule 163 descending truth table: ' + ''.join(map(str, table)))
    print('Fidelity: A266753 source terms n=0..24 and A392003 terms m=1..34 agree')
    print(f'Independent simulators: {cell_checks} cell comparisons, including both boundaries and exterior')
    print(f'Row invariant, encoding, closed form and B(2n)=A(n): n=1..{LIMIT}')
    print(f'First recurrence: n=2..{LIMIT}; second recurrence: n=3..{LIMIT}')
    print(f'Generating-function product: coefficients n=0..{LIMIT}')
    print('All arithmetic is integer or Fraction; finite checks are not a proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    run()
