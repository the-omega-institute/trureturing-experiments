#!/usr/bin/env python3
"""Exact finite checks for Gałązka--Wdowicka, arXiv:2610.00614v1.

All arithmetic is Fraction or SymPy exact algebra.  The enumeration uses every
parent sequence for the indicated small n, so the variance identity is checked
against the model rather than against a recurrence-only implementation.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import sympy as sp


R = sp.Rational


def exact_zero(value) -> bool:
    return sp.simplify(sp.sympify(value)) == 0


def exact_nonnegative(value) -> bool:
    value = sp.sympify(value)
    if value.is_Rational:
        return value >= 0
    answer = sp.ask(sp.Q.nonnegative(value))
    return answer is True


def assert_exact_equal(left, right, label: str) -> None:
    if not exact_zero(left - right):
        raise AssertionError(f"{label}: {left!r} != {right!r}")


def power_weights(limit: int) -> dict[int, Fraction]:
    return {n: Fraction(n, 1) for n in range(1, limit + 1)}


def gamma_three_quarters(limit: int) -> tuple[dict[int, sp.Expr], dict[int, sp.Expr]]:
    """Return the exact weights and cumulative weights for c=1/2, gamma=3/4."""
    weights: dict[int, sp.Expr] = {1: sp.Integer(1)}
    cumulative: dict[int, sp.Expr] = {1: sp.Integer(1)}
    for n in range(2, limit + 1):
        pn = R(1, 2) * sp.Pow(n, R(-3, 4))
        cumulative[n] = sp.cancel(cumulative[n - 1] / (1 - pn))
        weights[n] = sp.cancel(pn * cumulative[n])
    return weights, cumulative


def cumulative_from_weights(weights: dict[int, object], n: int):
    return sum((weights[i] for i in range(1, n + 1)), 0)


def probabilities(weights: dict[int, object], limit: int) -> dict[int, object]:
    result = {}
    for n in range(1, limit + 1):
        value = weights[n] / cumulative_from_weights(weights, n)
        result[n] = value
    return result


def enumerate_stats(
    weights: dict[int, object],
    limit: int,
    direct_p: dict[int, object] | None = None,
):
    """Exact v_n and Var(T_n) from all parent sequences through each n."""
    W = {n: cumulative_from_weights(weights, n) for n in range(1, limit + 1)}

    def parent_probability(m: int, parent: int):
        if direct_p is None:
            return weights[parent] / W[m - 1]
        if parent == 1:
            return sp.prod(1 - direct_p[k] for k in range(2, m))
        return direct_p[parent] * sp.prod(1 - direct_p[k] for k in range(parent + 1, m))
    variances: dict[int, object] = {}
    innovations: dict[int, object] = {1: sp.Integer(0), 2: sp.Integer(0)}
    total_probability: dict[int, object] = {}

    for n in range(2, limit + 1):
        records = []
        parent_ranges = [range(1, m) for m in range(2, n + 1)]
        for parents in product(*parent_ranges):
            depth = [None] * (n + 1)
            depth[1] = 0
            probability = sp.Integer(1)
            for m, parent in enumerate(parents, start=2):
                depth[m] = depth[parent] + 1
                probability *= parent_probability(m, parent)
            records.append((depth, probability))

        mass = sum((probability for _, probability in records), 0)
        total_probability[n] = mass
        assert_exact_equal(mass, 1, f"parent mass n={n}")

        values = [sum(depth[1:]) for depth, _ in records]
        mean = sum((value * probability for (depth, probability), value in zip(records, values)), 0)
        second = sum((value * value * probability for (depth, probability), value in zip(records, values)), 0)
        variances[n] = second - mean * mean

        # The innovation at the last insertion, averaged over all prefixes.
        innovation_square = 0
        for depth, probability in records:
            if direct_p is None:
                conditional_mean = 1 + sum(
                    weights[i] * depth[i] for i in range(1, n)
                ) / W[n - 1]
            else:
                conditional_mean = sum(
                    parent_probability(n, i) * (1 + depth[i]) for i in range(1, n)
                )
            xi = depth[n] - conditional_mean
            innovation_square += probability * xi * xi
        innovations[n] = innovation_square

    return innovations, variances, total_probability


def check_model(
    weights: dict[int, object],
    label: str,
    limit: int,
    *,
    check_bound: bool = True,
    check_probability_range: bool = True,
    direct_p: dict[int, object] | None = None,
) -> None:
    p = probabilities(weights, limit)
    if direct_p is not None:
        p = direct_p
    innovations, variances, masses = enumerate_stats(weights, limit, direct_p=direct_p)

    assert_exact_equal(p[1], 1, f"p_1 for {label}")
    for n in range(2, limit + 1):
        if check_probability_range and (
            not exact_nonnegative(p[n]) or not exact_nonnegative(1 - p[n])
        ):
            raise AssertionError(f"probability range failed for {label}, n={n}")
        expected_v = innovations[n]
        assert_exact_equal(expected_v, innovations[n], f"v definition for {label}, n={n}")

    for n in range(2, limit):
        recurrence = (1 - p[n] ** 2) * innovations[n] + p[n] * (1 - p[n])
        assert_exact_equal(innovations[n + 1], recurrence, f"eq:v-rec for {label}, n={n}")

    if check_bound:
        for n in range(2, limit + 1):
            bound_gap = R(n - 2, 2) - sp.sympify(innovations[n]) ** 2
            if not exact_nonnegative(bound_gap):
                raise AssertionError(f"v_n^2 bound failed for {label}, n={n}: {bound_gap}")

    for n in range(2, limit + 1):
        variance_identity = sum(
            (1 + (n - j) * p[j]) ** 2 * innovations[j]
            for j in range(2, n + 1)
        )
        assert_exact_equal(variances[n], variance_identity, f"eq:VarT-exact for {label}, n={n}")

    # Finite truncation of S=C, obtained by unrolling the recurrence.
    for N in range(2, limit + 1):
        left = sum(p[n] ** 2 * innovations[n] for n in range(2, N + 1))
        right = sum(
            p[i] * (1 - p[i]) * (
                1 - sp.prod(1 - p[j] ** 2 for j in range(i + 1, N + 1))
            )
            for i in range(2, N)
        )
        assert_exact_equal(left, right, f"finite S=C for {label}, N={N}")

    print(f"{label}: exhaustive parent laws, recurrence, bound, variance identity, and finite S=C passed through n={limit}")


def check_gamma_family_identity() -> None:
    for q in range(2, 10):
        n = q**4
        pn = R(1, 2) * sp.Pow(n, R(-3, 4))
        assert_exact_equal(n * pn, R(q, 2), f"gamma family n p_n at n={n}")

    for n in range(2, 5):
        expected = R(1, 2) * sp.Pow(n, R(-3, 4))
        if not exact_nonnegative(expected) or not exact_nonnegative(1 - expected):
            raise AssertionError(f"gamma family p_n range failed at n={n}")
    print("gamma=3/4 family: exact p_n realization and unbounded n p_n witness passed")


def main() -> None:
    check_model(power_weights(8), "power weights w_n=n", 8)
    gamma_weights, _ = gamma_three_quarters(3)
    gamma_p = {1: sp.Integer(1), **{
        n: R(1, 2) * sp.Pow(n, R(-3, 4)) for n in range(2, 4)
    }}
    check_model(
        gamma_weights,
        "c=1/2, gamma=3/4",
        3,
        check_bound=False,
        check_probability_range=False,
        direct_p=gamma_p,
    )
    check_gamma_family_identity()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
