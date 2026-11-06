#!/usr/bin/env python3
"""Exact interval check for Manjunath, arXiv:2410.15426, Section 7.

The intervals are one-dimensional polyhedra.  Multiplication is Minkowski
sum, so a power scales both endpoints; addition is the convex hull of the
union.  All arithmetic is over Fraction.
"""

from fractions import Fraction
from itertools import combinations


Interval = tuple[Fraction, Fraction]


def add(a: Interval, b: Interval) -> Interval:
    return (a[0] + b[0], a[1] + b[1])


def scale(a: Interval, n: int) -> Interval:
    return (n * a[0], n * a[1])


def hull(intervals: list[Interval]) -> Interval:
    return (min(a for a, _ in intervals), max(b for _, b in intervals))


def product_coefficients(factors: list[Interval]) -> list[Interval]:
    degree = len(factors)
    coefficients: list[Interval] = []
    for i in range(degree + 1):
        choices = [
            (sum((factors[j][0] for j in subset), Fraction(0)),
             sum((factors[j][1] for j in subset), Fraction(0)))
            for subset in combinations(range(degree), degree - i)
        ]
        coefficients.append(hull(choices) if choices else (Fraction(0), Fraction(0)))
    return coefficients


# The ordinary quartic discriminant, in variables a,b,c,d,e, with signs
# suppressed by Newton-polytope polyhedralisation.  The 16 exponent vectors
# are the complete support of the discriminant.
DISCRIMINANT_EXPONENTS = [
    (3, 0, 0, 0, 3), (2, 1, 0, 1, 2), (2, 0, 2, 0, 2),
    (2, 0, 1, 2, 1), (2, 0, 0, 4, 0), (1, 2, 1, 0, 2),
    (1, 2, 0, 2, 1), (1, 1, 2, 1, 1), (1, 1, 1, 3, 0),
    (1, 0, 4, 0, 1), (1, 0, 3, 2, 0), (0, 4, 0, 0, 2),
    (0, 3, 1, 1, 1), (0, 3, 0, 3, 0), (0, 2, 3, 0, 1),
    (0, 2, 2, 2, 0),
]


def monomial_interval(coefficients: list[Interval], exponents: tuple[int, ...]) -> Interval:
    result = (Fraction(0), Fraction(0))
    for coefficient, exponent in zip(coefficients, exponents):
        result = add(result, scale(coefficient, exponent))
    return result


def discriminant_evaluation(coefficients: list[Interval]) -> tuple[Interval, list[Interval]]:
    terms = [monomial_interval(coefficients, exponents) for exponents in DISCRIMINANT_EXPONENTS]
    return hull(terms), terms


def slope_at_breakpoints(breakpoints: list[Fraction], slopes: list[int]) -> dict[Fraction, tuple[int, int]]:
    return {point: (slopes[index], slopes[index + 1]) for index, point in enumerate(breakpoints)}


def main() -> None:
    factors = [(Fraction(2), Fraction(3)), (Fraction(2), Fraction(4)),
               (Fraction(-1), Fraction(1)), (Fraction(0), Fraction(1))]
    q = product_coefficients(factors)
    assert q == [(Fraction(3), Fraction(9)), (Fraction(1), Fraction(8)),
                 (Fraction(-1), Fraction(7)), (Fraction(-1), Fraction(4)),
                 (Fraction(0), Fraction(0))]
    # a,b,c,d,e = Q4,Q3,Q2,Q1,Q0.
    evaluation, terms = discriminant_evaluation([q[4], q[3], q[2], q[1], q[0]])
    assert evaluation == (Fraction(-2), Fraction(38))
    endpoint_hits = [sum(term[0] == evaluation[0] for term in terms),
                     sum(term[1] == evaluation[1] for term in terms)]
    assert endpoint_hits == [2, 2]
    assert terms[-2] == evaluation and terms[-1] == evaluation
    print("fidelity: Q0..Q4 =", q)
    print("fidelity: quartic discriminant support rows =", len(terms))
    print("root evaluation:", evaluation, "endpoint multiplicities:", endpoint_hits)

    # Affine covariance of the same family.  For P'_j=t+lambda P_j,
    # Q'_i=(4-i)t+lambda Q_i, and the discriminant interval is predicted to be
    # [12t-2lambda, 12t+38lambda].
    for t, lam in [(Fraction(-3), Fraction(2)), (Fraction(0), Fraction(1)),
                   (Fraction(5), Fraction(3)), (Fraction(7, 2), Fraction(5, 2))]:
        q_affine = [add(((4 - i) * t, (4 - i) * t),
                         (lam * q[i][0], lam * q[i][1])) for i in range(5)]
        affine_eval, _ = discriminant_evaluation(
            [q_affine[4], q_affine[3], q_affine[2], q_affine[1], q_affine[0]])
        assert affine_eval == (12 * t - 2 * lam, 12 * t + 38 * lam)
    print("affine family: checked four rational (t, lambda) instances")

    # Endpoint slope obstruction.  A lower endpoint with three distinct
    # integer active exponents needs a slope drop of at least two; the only
    # such breakpoint of F is u=2.  The upper endpoint analog for G has a
    # slope rise of at least two only at v=1.
    f_jumps = slope_at_breakpoints([Fraction(-1), Fraction(0), Fraction(2)],
                                    [4, 3, 2, 0])
    g_jumps = slope_at_breakpoints([Fraction(1), Fraction(3), Fraction(4)],
                                    [0, 2, 3, 4])
    lower_candidates = [x for x, (left, right) in f_jumps.items() if left - right >= 2]
    upper_candidates = [x for x, (left, right) in g_jumps.items() if right - left >= 2]
    assert lower_candidates == [Fraction(2)]
    assert upper_candidates == [Fraction(1)]
    assert lower_candidates[0] > upper_candidates[0]
    print("slope obstruction: F forces u=2 and G forces v=1")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
