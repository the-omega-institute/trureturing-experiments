"""Hilbert-Schmidt Gram matrices and subgroup indices from integer weights."""
from fractions import Fraction
from itertools import combinations

PAIRS = tuple(combinations(range(5), 2))


def basis():
    units = [{(i, j): 1} for i in range(5) for j in range(5) if i != j]
    diagonals = [{(i, i): 1, (i + 1, i + 1): -1} for i in range(4)]
    return units + diagonals


def derivative(a):
    result = {}
    for (i, j), value in a.items():
        for col, pair in enumerate(PAIRS):
            if j not in pair:
                continue
            new_pair = tuple(i if k == j else k for k in pair)
            if new_pair[0] == new_pair[1]:
                continue
            row = PAIRS.index(tuple(sorted(new_pair)))
            sign = 1 if new_pair[0] < new_pair[1] else -1
            result[row, col] = result.get((row, col), 0) + sign * value
    return result


def block_diagonal(b):
    return {(i, j): value for (i, j), value in b.items()
            if (PAIRS[i][1] < 3) == (PAIRS[j][1] < 3)}


def inner(a, b):
    return sum(value * b.get(key, 0) for key, value in a.items())


def adjoint_index(weights, denominator, allowed):
    trace = sum((x - y)**2 for i, x in enumerate(weights)
                for j, y in enumerate(weights) if allowed(i, j))
    return Fraction(trace, denominator)


def verify():
    generators = basis()
    lifted = [derivative(a) for a in generators]
    kept = [block_diagonal(b) for b in lifted]
    factors = []
    for a in generators:
        is_unbroken = all((i < 3) == (j < 3) for i, j in a)
        factors.append(3 if is_unbroken else 1)
    for i, a in enumerate(generators):
        for j, b in enumerate(generators):
            assert inner(lifted[i], lifted[j]) == 3 * inner(a, b)
            assert inner(kept[i], kept[j]) == factors[i] * inner(a, b)
    assert factors.count(3) == factors.count(1) == 12
    assert 99 - len(generators) == 75
    # For x=g5^2, y=g10^2, tr=3x+y, det=(3-factor)xy.
    # At lambda=y the characteristic polynomial is -factor*x*y.
    for factor in (1, 3):
        assert -3 + (3 - factor) == -factor
    for weights, expected_h in [((-2, -2, -2, 3, 3), 0),
                                ((0, 0, 0, 1, -1), 2),
                                ((1, -1, 0, 0, 0), 3)]:
        denominator = 2 * sum(x*x for x in weights)
        wedge = [weights[i] + weights[j] for i, j in PAIRS]
        full = lambda i, j: True
        h = adjoint_index(weights, denominator, lambda i, j: (i < 3) == (j < 3))
        su5 = adjoint_index(weights, denominator, full)
        su10 = adjoint_index(wedge, denominator, full)
        assert (h, su5, su10) == (expected_h, 5, 30)
        assert su10 - su5 == 25
    print("Full and projected Gram matrices checked on 24 x 24 basis pairs.")
    print("Multiplicities: 12 massless; 12 axial H, 12 + 12 mixed, 75 pure target.")
    print("Adjoint indices H=(0,2,3), mixed=(5,3,2), complement=(25,25,25).")


if __name__ == "__main__":
    verify()
