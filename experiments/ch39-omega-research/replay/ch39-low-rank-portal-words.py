"""Exact compressed spectral weights and the first distinguishing source word."""
from fractions import Fraction
from math import prod

from sympy import Matrix, eye, sqrt, zeros


def compressed_moments(n, rotated):
    nodes = range(1, n + 1)
    weights = []
    identity = eye(2)
    for a in nodes:
        derivative = prod(a - b for b in nodes if a != b)
        if rotated:
            direction = Matrix([[a*a - 1, 2*a], [2*a, 1 - a*a]])
        else:
            direction = Matrix.diag(a*a + 1, -a*a - 1)
        weight = Fraction(a*a + 1, abs(derivative))*identity + direction/derivative
        assert weight == weight.T and weight.det() == 0 and weight.trace() > 0
        weights.append(weight)
    total = sum(weights, zeros(2))
    assert total == total[0, 0]*identity
    scale = 1/total[0, 0]
    return scale, [
        scale*sum((a**k * w for a, w in zip(nodes, weights)), zeros(2))
        for k in range(2*n - 3)
    ]


def check_words(n):
    order = 2*n - 4
    scale, first = compressed_moments(n, False)
    other_scale, second = compressed_moments(n, True)
    assert scale == other_scale
    # Every word with projections expands in I and P. When at least one P
    # occurs, its cyclic trace is the trace of a composition of positive gaps.
    count = 0

    def visit(total, left, right):
        nonlocal count
        if total >= 3:
            assert left.trace() == right.trace()
            count += 1
        for gap in range(1, order - total):
            visit(total + gap, left*first[gap], right*second[gap])

    visit(0, eye(2), eye(2))
    signal = (second[n - 2]**2 - first[n - 2]**2).trace()
    assert signal == 8*scale**2
    print(f"N={n}: {count} subcritical compositions; first witness at {order}, difference {signal}")


if __name__ == "__main__":
    first, second = Matrix.diag(1, 0), eye(2)/sqrt(2)  # two choices of BB*
    assert (first**2).trace() == (second**2).trace() == 1
    assert (first**3).trace() == 1 and (second**3).trace() == 1/sqrt(2)
    for dimension in (6, 7, 8):
        check_words(dimension)
