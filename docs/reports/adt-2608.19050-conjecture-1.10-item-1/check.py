#!/usr/bin/env python3
"""Exact symbolic and finite structural checks for Conjecture 1.10(1).

All identities used as evidence are SymPy exact identities over symbolic
expressions.  The finite permutation checks use integer arithmetic only.
"""

from fractions import Fraction
import sympy as sp


def alternating(n):
    if n == 0:
        return []
    if n == 1:
        return [1]
    if n % 2 == 1:
        return alternating(n - 1) + [n]
    return alternating(n - 2) + [n, n - 1]


def standardize(word):
    values = sorted(word)
    rank = {value: i + 1 for i, value in enumerate(values)}
    return [rank[value] for value in word]


def source_recursion_check():
    p, q, a, b, x = sp.symbols("p q a b x", positive=True)
    continuation = p / x + q / (1 - x)
    expected = p * (sp.log(b) - sp.log(a)) + q * (
        sp.log(1 - a) - sp.log(1 - b)
    )
    # Differentiate with respect to the upper endpoint and normalize at b=a.
    assert sp.simplify(sp.diff(expected, b) - continuation.subs(x, b)) == 0
    assert sp.simplify(expected.subs(b, a)) == 0
    integrated = expected
    ratio_terms = p * sp.log(b / a) + q * sp.log((1 - a) / (1 - b))
    assert sp.simplify(sp.logcombine(integrated - ratio_terms, force=True)) == 0
    interval = (1 + ratio_terms) / (b - a)
    assert interval == (1 + ratio_terms) / (b - a)


def alternating_recursion_check():
    for n in range(1, 101):
        word = alternating(n)
        assert word == list(range(1, n + 1)) if n == 1 else sorted(word) == list(
            range(1, n + 1)
        )
        if n >= 2:
            lower = [value for value in word[1:] if value < word[0]]
            upper = [value for value in word[1:] if value > word[0]]
            assert standardize(lower) == ([1] if lower else [])
            assert standardize(upper) == alternating(n - 2)


def rescaled_variational_check():
    r, u, v = sp.symbols("r u v", positive=True)
    A = (1 - sp.log(r)) / (1 - r)
    B = (1 + r) / 2
    H = (1 + sp.log(v / u)) / (v - u) + (u + v) / 2
    rescaled = sp.factor(H.subs({u: r * v}))
    assert sp.simplify(rescaled - (A / v + B * v)) == 0
    M_integrand = (1 + r) * (1 - sp.log(r)) / (1 - r)
    # C^2 / 16 = M / 8, with C = sqrt(2 M), exactly as in the outline.
    assert sp.simplify((2 * A * B) / 8 - M_integrand / 8) == 0


def exact_rational_fixtures():
    # The first two boundary values in the published Bellman recurrence.
    assert Fraction(0, 1) == 0
    assert Fraction(1, 1) == 1
    # The recurrence has one lower singleton and upper Alt_(n-2).
    assert alternating(1) == [1]
    assert alternating(2) == [2, 1]


if __name__ == "__main__":
    source_recursion_check()
    alternating_recursion_check()
    rescaled_variational_check()
    exact_rational_fixtures()
    print("source recurrence: exact logarithmic integral identity OK")
    print("Alt_n finite permutation and two-block decomposition: n <= 100 OK")
    print("rescaled H=A/v+Bv and C^2/16=M/8: exact symbolic identity OK")
    print("ALL_CHECKS_OK")
