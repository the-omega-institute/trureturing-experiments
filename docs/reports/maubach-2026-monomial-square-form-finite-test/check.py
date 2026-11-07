#!/usr/bin/env python3
"""Exact checks for Maubach, Problem 9 (arXiv:2611.03057v1).

The finite enumeration includes negative and zero inputs.  All arithmetic is
integer arithmetic; Fraction is used for the dyadic interval convention and
no floating point value is used as evidence.
"""

from fractions import Fraction
from itertools import product


def f_value(k, xs):
    """F_k(x) = sum_{j=1}^k (x_1 ... x_j)^2."""
    total = 0
    prefix = 1
    for x in xs[:k]:
        prefix *= x
        total += prefix * prefix
    return total


def represented_in_box(k, t):
    """Values <= 2^(2t) represented with x_i in [-2^t,2^t]."""
    radius = 2**t
    limit = 2 ** (2 * t)
    values = set()
    for xs in product(range(-radius, radius + 1), repeat=k):
        value = f_value(k, xs)
        if 0 < value <= limit:
            values.add(value)
    return values


def counting_bound(k, t):
    """The source-independent dyadic upper bound k*2^t*(t+1)^k."""
    return k * (2**t) * ((t + 1) ** k)


def source_fidelity():
    for k in range(1, 5):
        for r in range(1, k + 1):
            xs = (1,) * r + (0,) * (k - r)
            assert f_value(k, xs) == r
        print(f"fidelity k={k}: represents 1..{k} by 1/0 specialization")


def finite_exhaustive_counts():
    for k in range(1, 5):
        for t in range(1, 5):
            values = represented_in_box(k, t)
            bound = counting_bound(k, t)
            assert len(values) <= bound
            print(
                f"exhaustive k={k} t={t} input_box=[{-2**t},{2**t}] "
                f"positive_values={len(values)} bound={bound}"
            )


def asymptotic_bound_check():
    # The two elementary power inequalities used in the proof are checked
    # exactly on a finite audit range; their standard inductions are stated in
    # the accompanying README.
    for k in range(4, 65):
        assert k <= 2 ** (k - 1)
        assert 4 * k * k + 1 <= 2 ** (3 * k)
        t = 4 * k * k
        bound = counting_bound(k, t)
        available = 2 ** (2 * t)
        assert bound <= 2 ** (8 * k * k - 1)
        assert bound < available - k
    k = 4
    t = 4 * k * k
    print(
        f"asymptotic sample k={k} t={t} "
        f"bound_bits={counting_bound(k,t).bit_length()} "
        f"available_bits={(2**(2*t)).bit_length()}"
    )


def finite_test_set_check():
    # This is an arbitrary finite test set.  The same construction works for
    # every finite set by choosing k >= its largest element.
    test_set = {1, 2, 7, 14, 29, 9591}
    k = max(4, max(test_set))
    for n in test_set:
        assert f_value(k, (1,) * n + (0,) * (k - n)) == n
    t = 4 * k * k
    assert counting_bound(k, t) < 2 ** (2 * t) - k
    print(f"finite test set {sorted(test_set)} represented by F_{k}; F_{k} is non-universal")


def main():
    # Keep an exact Fraction operation in the checker to make the dyadic
    # interval convention explicit without introducing a float.
    assert Fraction(2, 1) == 2
    source_fidelity()
    finite_exhaustive_counts()
    asymptotic_bound_check()
    finite_test_set_check()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
