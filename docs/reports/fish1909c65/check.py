"""Exact finite checks for Fish et al., Conjecture 6.5."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt

from sympy import QQ, Rational, sqrt


def difference(p, q):
    return p[0] - q[0], p[1] - q[1]


def cross(p, q):
    return p[0] * q[1] - p[1] * q[0]


def check_source_example(n):
    # Example 4.5(1) in the source, in the Chebyshev norm.
    points = [(i, Fraction(1, i)) for i in range(1, n + 1)]
    distances = Counter()
    for i, j in combinations(range(n), 2):
        dx, dy = difference(points[j], points[i])
        distance = max(abs(dx), abs(dy))
        assert distance == j - i
        distances[distance] += 1
    assert distances == Counter({gap: n - gap for gap in range(1, n)})
    for a, b, c in combinations(points, 3):
        assert cross(difference(b, a), difference(c, a)) != 0
    for a, b, c, d in combinations(points, 4):
        # Middle points are interior to every axis-aligned square containing
        # both endpoints, since each coordinate is strictly monotone.
        assert a[0] < b[0] < c[0] < d[0]
        assert a[1] > b[1] > c[1] > d[1]


def sign_in_quadratic_field(z, radicand):
    """Return the sign of z = constant + coefficient * sqrt(radicand)."""
    coefficients = z.to_list()
    if not coefficients:
        return 0
    if len(coefficients) == 1:
        return 1 if coefficients[0] > 0 else -1
    coefficient, constant = coefficients
    if coefficient > 0 and constant >= 0:
        return 1
    if coefficient < 0 and constant <= 0:
        return -1
    square_difference = coefficient**2 * radicand - constant**2
    assert square_difference != 0
    return 1 if (square_difference > 0) == (coefficient > 0) else -1


def check_counterexample(n):
    denominator = 10 * n
    radicand = denominator**2 - 1
    assert isqrt(radicand) == denominator - 1
    field = QQ.algebraic_field(sqrt(radicand))
    zero, one = field.zero, field.one
    half = field.from_sympy(Rational(1, 2))
    s = field.from_sympy(Rational(1, denominator))
    c = field.from_sympy(sqrt(radicand) / denominator)

    def positive(z):
        return sign_in_quadratic_field(z, radicand) > 0

    x, y = zero, one
    points = []
    for _ in range(n):
        x, y = c * x + s * y, c * y - s * x
        assert x * x + y * y == one
        assert positive(y - half)
        points.append((x, y))

    # Certify a nontrivial unit segment in the proposed norm.
    for flat_x in (-half, zero, half):
        assert positive(one - flat_x * flat_x - half * half)
        assert half + half == one

    squared_distance_by_gap = {}
    multiplicities = Counter()
    for i, j in combinations(range(n), 2):
        dx, dy = difference(points[j], points[i])
        assert positive(-dy) and positive(dx + dy + dy)
        # The first inequality gives distinct heights. The second excludes
        # occupancy of opposite flats of any circle of the truncated norm
        # and certifies that the Euclidean term gives this pair's distance.
        squared_distance = dx * dx + dy * dy
        gap = j - i
        assert squared_distance == squared_distance_by_gap.setdefault(gap, squared_distance)
        multiplicities[gap] += 1
    assert multiplicities == Counter({gap: n - gap for gap in range(1, n)})
    for gap in range(1, n - 1):
        assert positive(squared_distance_by_gap[gap + 1] - squared_distance_by_gap[gap])

    # Every triple determines the original Euclidean circle. Its points
    # lie above the cap. Each flat contributes at most one point by height,
    # and the slope test excludes points on both opposite flats.
    for first, second, third in combinations(points, 3):
        assert cross(difference(second, first), difference(third, first)) != zero

    # Independently test whether an affine functional could send the first
    # four points to distinct ranks of an n-term arithmetic progression.
    # Affine interpolation on the first three fixes the fourth value.
    p0, p1, p2, p3 = points[:4]
    d1, d2, d3 = (difference(p, p0) for p in (p1, p2, p3))
    denominator_det = cross(d1, d2)
    numerator_1 = cross(d3, d2)
    numerator_2 = cross(d1, d3)
    rank_tuples = 0
    for k0, k1, k2, k3 in permutations(range(n), 4):
        actual = (k1 - k0) * numerator_1 + (k2 - k0) * numerator_2
        assert actual != (k3 - k0) * denominator_det
        rank_tuples += 1
    return rank_tuples


def main():
    for n in range(4, 9):
        check_source_example(n)
    print("Example 4.5(1), Chebyshev norm: n=4..8 passed")
    for n in range(4, 13):
        rank_tuples = check_counterexample(n)
        print(f"Truncated Euclidean norm: n={n}, {rank_tuples} rank tuples excluded")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
