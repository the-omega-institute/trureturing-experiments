"""Exact coloring checks for Clark–Elhamdadi–Saito–Yeatman Conjecture 5.2."""

from collections import Counter
from itertools import permutations, product
from math import prod


def check_trefoil():
    # Model C[4,1] as F_2[t]/(t^2+t+1), with bitwise XOR for addition.
    def multiply_by_t(x):
        return (x << 1) ^ (0b111 if x & 0b10 else 0)

    def crossing(pair):
        x, y = pair
        # Figure 1, right strand over: (x, y) maps to (y, x * y).
        return y, multiply_by_t(x ^ y) ^ y

    colors = [
        (x, y)
        for x, y in product(range(4), repeat=2)
        if crossing(crossing(crossing((x, y)))) == (x, y)
    ]
    trivial = sum(x == y for x, y in colors)
    assert (len(colors), trivial) == (16, 4)
    print("Source fidelity (trefoil, C[4,1]): 16 total = 4 trivial + 12 nontrivial")


def transpose(matrix):
    return list(map(list, zip(*matrix, strict=True)))


def determinant(matrix):
    size = len(matrix)
    result = 0
    for order in permutations(range(size)):
        inversions = sum(order[i] > order[j] for i in range(size) for j in range(i + 1, size))
        result += (-1) ** inversions * prod(matrix[i][order[i]] for i in range(size))
    return result


def seifert_matrix(p):
    upper = [[2, 1, 0], [0, 2, p], [0, 0, 2]]
    lower = [[entry - (i == j) for j, entry in enumerate(row)] for i, row in enumerate(transpose(upper))]
    matrix = [[0, 0, 0] + row for row in upper] + [row + [0, 0, 0] for row in lower]
    skew = [
        [x - y for x, y in zip(row, transposed, strict=True)]
        for row, transposed in zip(matrix, transpose(matrix), strict=True)
    ]
    assert determinant(skew) == 1
    return matrix


def action_parameters(p):
    modulus = p * p
    return (pow(2, -1, modulus), pow(4, -1, modulus)), (2, -1)


def rank_with_unit_pivots(matrix, modulus, prime):
    """Eliminate units by invertible operations and return the trailing block."""
    block = [[entry % modulus for entry in row] for row in matrix]
    pivots = 0
    while block and block[0]:
        pivot = next(
            ((i, j) for i, row in enumerate(block) for j, entry in enumerate(row)
             if entry % prime != 0),
            None,
        )
        if pivot is None:
            break
        i, j = pivot
        block[0], block[i] = block[i], block[0]
        for row in block:
            row[0], row[j] = row[j], row[0]
        scale = pow(block[0][0], -1, modulus)
        block = [
            [(entry - row[0] * scale * head) % modulus
             for entry, head in zip(row[1:], block[0][1:], strict=True)]
            for row in block[1:]
        ]
        pivots += 1
    return pivots, block


def mixed_kernel_size(left, right, p, parameters):
    """Count the kernel of U left − right using a lift to Z/p² coordinates."""
    u, off = parameters
    lifted = []
    for left_row, right_row in zip(left, right, strict=True):
        long_row, short_row = [], []
        for a, b in zip(left_row, right_row, strict=True):
            diagonal = u * a - b
            long_row.extend((diagonal, p * off * a))
            short_row.extend((0, p * diagonal))
        lifted.extend((long_row, short_row))

    pivots, residual = rank_with_unit_pivots(lifted, p * p, p)
    rank, trailing = rank_with_unit_pivots(
        [[entry // p for entry in row] for row in residual], p, p
    )
    assert all(entry == 0 for row in trailing for entry in row)
    size = p ** (2 * (len(lifted[0]) - pivots) - rank)
    # Each original short coordinate has p lifts.
    assert size % p ** len(left[0]) == 0
    return size // p ** len(left[0])


def enumerated_kernel_size(left, right, p, parameters):
    """Count original M³ inputs by matching the third column contribution."""
    u, off = parameters
    moduli = (p * p, p) * len(left)
    points = list(product(range(p * p), range(p)))
    columns = []
    for a_col, b_col in zip(transpose(left), transpose(right), strict=True):
        columns.append([
            tuple(
                coordinate
                for a, b in zip(a_col, b_col, strict=True)
                for coordinate in (((u * a - b) * x + p * off * a * y) % (p * p),
                                   (u * a - b) * y % p)
            )
            for x, y in points
        ])
    first, second, third = columns
    multiplicities = Counter(third)
    return sum(
        multiplicities[tuple((-x - y) % modulus for x, y, modulus in zip(a, b, moduli, strict=True))]
        for a, b in product(first, second)
    )


def check_target(p, parameters):
    points = set(product(range(p * p), range(p)))

    def evaluate(action, x, y):
        u, off = action
        return (u * x + p * off * y) % (p * p), u * y % p

    forward, inverse = parameters
    assert all(evaluate(inverse, *evaluate(forward, x, y)) == (x, y) for x, y in points)
    displacement = [(x - a, y - b) for x, y in points for a, b in [evaluate(forward, x, y)]]
    assert {(a % (p * p), b % p) for a, b in displacement} == points


def check_family():
    for p in (5, 7, 11, 13, 17, 19, 31):
        matrix = seifert_matrix(p)
        matrix_t = transpose(matrix)
        blocks = [
            ([row[3:] for row in matrix[:3]], [row[3:] for row in matrix_t[:3]]),
            ([row[:3] for row in matrix[3:]], [row[:3] for row in matrix_t[3:]]),
        ]
        parameters = action_parameters(p)
        check_target(p, parameters)
        kernels = [tuple(mixed_kernel_size(a, b, p, action) for a, b in blocks) for action in parameters]
        assert kernels == [(p**4, 1), (1, p**5)]
        if p == 5:
            enumerated = [
                tuple(enumerated_kernel_size(a, b, p, action) for a, b in blocks)
                for action in parameters
            ]
            assert enumerated == kernels
            print(f"p=5 original-block enumeration ({p**9:,} inputs per block): {enumerated}")
        totals = tuple(p**3 * prod(pair) for pair in kernels)
        assert totals == (p**7, p**8)
        print(f"p={p}: kernels {kernels}, coloring totals {totals}")


if __name__ == "__main__":
    check_trefoil()
    check_family()
    print("ALL_CHECKS_OK")
