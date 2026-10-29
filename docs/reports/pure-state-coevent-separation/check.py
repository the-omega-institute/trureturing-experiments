#!/usr/bin/env python3
"""Exact finite checks for Wallden's coevent-separation appendix.

The script uses only integers, Fraction, and SymPy exact expressions.  It
checks the eight branch amplitudes, the minimal null events, and the primitive
support classifications for the real protocol in the source and for the
phase-regularized protocol nominated for the settlement.  Extra zero branch
atoms model the orthogonal complement used in the finite-dimensional lift.
"""

from fractions import Fraction
from itertools import product

import sympy as sp


I = sp.I
SQRT2 = sp.sqrt(2)


def vec(block, scalar):
    return (block, scalar)


def exact_zero(x):
    return x == (Fraction(0), Fraction(0))


def add_scalar(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg_scalar(x):
    return (-x[0], -x[1])


def event_zero(branches, mask):
    totals = [(Fraction(0), Fraction(0)), (Fraction(0), Fraction(0))]
    for i, (block, scalar) in enumerate(branches):
        if mask & (1 << i):
            totals[block] = add_scalar(totals[block], scalar)
    return all(exact_zero(x) for x in totals)


def null_masks(branches):
    return {
        mask
        for mask in range(1 << len(branches))
        if event_zero(branches, mask)
    }


def minimal_nonempty(nulls):
    answer = set()
    for mask in nulls:
        if mask == 0:
            continue
        proper_nonempty_null = any(
            other != 0
            and other != mask
            and (other & ~mask) == 0
            for other in nulls
        )
        if not proper_nonempty_null:
            answer.add(mask)
    return answer


def primitive_supports(branches):
    nulls = null_masks(branches)

    def preclusive(mask):
        return not any((mask & ~zero) == 0 for zero in nulls)

    result = set()
    for mask in range(1, 1 << len(branches)):
        if not preclusive(mask):
            continue
        has_preclusive_proper_subset = any(
            other != 0
            and other != mask
            and (other & ~mask) == 0
            and preclusive(other)
            for other in range(1, 1 << len(branches))
        )
        if not has_preclusive_proper_subset:
            result.add(mask)
    return result


def source_branches(c, s):
    u = ((c - s) / 2, Fraction(0))
    v = ((c + s) / 2, Fraction(0))
    return [
        vec(0, u), vec(0, u), vec(0, v), vec(0, neg_scalar(v)),
        vec(1, u), vec(1, neg_scalar(u)), vec(1, v), vec(1, v),
    ]


def phase_branches(c, s):
    u = (c / 2, -s / 2)
    v = (c / 2, s / 2)
    return [
        vec(0, u), vec(0, u), vec(0, v), vec(0, neg_scalar(v)),
        vec(1, u), vec(1, neg_scalar(u)), vec(1, v), vec(1, v),
    ]


def state_two_branches():
    a = (Fraction(1, 2), Fraction(0))
    return [
        vec(0, a), vec(0, a), vec(0, a), vec(0, neg_scalar(a)),
        vec(1, a), vec(1, neg_scalar(a)), vec(1, a), vec(1, a),
    ]


def assert_classification(branches, expected_min_null, expected_primitive):
    nulls = null_masks(branches)
    assert minimal_nonempty(nulls) == expected_min_null, (minimal_nonempty(nulls), expected_min_null)
    assert primitive_supports(branches) == expected_primitive, (
        primitive_supports(branches), expected_primitive
    )


def check_protocol(c, s, phase):
    branches = phase_branches(c, s) if phase else source_branches(c, s)
    expected_min_null = {0b00001100, 0b00110000}
    expected_c1 = {1, 2, 64, 128}
    assert_classification(branches, expected_min_null, expected_c1)

    branches_two = state_two_branches()
    expected_two_min_null = {9, 10, 12, 48, 96, 160}
    expected_c2 = {3, 5, 6, 80, 144, 192}
    assert_classification(branches_two, expected_two_min_null, expected_c2)
    assert expected_c1.isdisjoint(expected_c2)

    # The full event list includes unions of minimal null events.  This checks
    # the exact event definition rather than only the minimal list quoted in
    # the paper.
    assert len(null_masks(branches_two)) == 16

    # Adding zero branch atoms is the finite analogue of adjoining I-P_S.
    for zero_count in (1, 2):
        for base in (branches, branches_two):
            extended = base + [vec(0, (Fraction(0), Fraction(0)))] * zero_count
            supports = primitive_supports(extended)
            assert all(mask < (1 << 8) for mask in supports)


def main():
    # Exact Pythagorean representatives, all away from the source's
    # tan(theta)=+/-1/3 exceptional values.
    pairs = [
        (Fraction(3, 5), Fraction(4, 5)),
        (Fraction(5, 13), Fraction(12, 13)),
        (Fraction(8, 17), Fraction(15, 17)),
    ]
    for c, s in pairs:
        c = Fraction(c)
        s = Fraction(s)
        assert sp.simplify(sp.Rational(c.numerator, c.denominator) ** 2 +
                           sp.Rational(s.numerator, s.denominator) ** 2 - 1) == 0
        check_protocol(c, s, phase=False)
        check_protocol(c, s, phase=True)

        # The phase construction removes all bounded integer cancellations:
        # r u + t v = 0 has only r=t=0 for the tested exact instances.
        u = (c / 2, -s / 2)
        v = (c / 2, s / 2)
        for r, t in product(range(-12, 13), repeat=2):
            if (r, t) != (0, 0):
                assert not exact_zero((r * u[0] + t * v[0], r * u[1] + t * v[1]))

    print("checked three exact normalized overlaps")
    print("checked source and phase protocols: null events, primitive supports, disjointness")
    print("checked one- and two-atom zero-vector extensions")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
