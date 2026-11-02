#!/usr/bin/env python3
"""Exact small checks for Cordoba-Caycedo--Garcia Conjecture 4.30.

All cardinality calculations use integers or Fraction.  The script checks the
primary decomposition formulas used in the paper and the bounded-width
p-torsion obstruction from the preregistered proof route.
"""
from __future__ import annotations

from fractions import Fraction
from math import prod
from typing import Iterable

Factor = tuple[int, int, int]  # (prime, exponent, multiplicity)


def vp(n: int, p: int) -> int:
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def group_order(factors: Iterable[Factor]) -> int:
    return prod(p ** (j * t) for p, j, t in factors)


def p_width(factors: Iterable[Factor], p: int) -> int:
    return max((t for q, _j, t in factors if q == p), default=0)


def p_rank(factors: Iterable[Factor], p: int) -> int:
    return sum(t for q, _j, t in factors if q == p)


def kernel_size(n: int, factors: Iterable[Factor]) -> int:
    """Cardinality of G[n] from the cyclic primary decomposition."""
    return prod(p ** (t * min(vp(n, p), j)) for p, j, t in factors)


def quotient_size(n: int, factors: Iterable[Factor]) -> int:
    """Cardinality of G/nG, equal to |G[n]| by the first isomorphism theorem."""
    return prod(p ** (t * min(vp(n, p), j)) for p, j, t in factors)


def direct_sum(*blocks: Factor) -> tuple[Factor, ...]:
    return tuple(blocks)


def check_decomposition_counts() -> None:
    g = direct_sum((2, 1, 2), (2, 2, 1), (3, 1, 1))
    assert group_order(g) == 2**4 * 3
    assert p_width(g, 2) == 2
    assert p_rank(g, 2) == 3
    assert kernel_size(12, g) == 2**4 * 3
    assert quotient_size(12, g) == kernel_size(12, g)

    # The two positive equations 2x=0 and 3x=0 have intersection G[gcd(2,3)]=G[1].
    assert kernel_size(2, g) == 2**3
    assert kernel_size(3, g) == 3
    assert kernel_size(1, g) == 1


def check_bounded_width_obstruction() -> None:
    # G_r = direct sum_{j=1}^r Z/2^j.  Its 2-width is one, while its total
    # 2-rank and 2-torsion both grow.  For N=2, the prescribed N-dimensional
    # scale is |G_r|^(1/2), and the exact ratio tends to zero.
    ratios: list[Fraction] = []
    for r in range(1, 13):
        g = tuple((2, j, 1) for j in range(1, r + 1))
        assert p_width(g, 2) == 1
        assert p_rank(g, 2) == r
        torsion = kernel_size(2, g)
        order = group_order(g)
        # |G_r|^(1/2) is an integer exactly when the exponent is even; use
        # the squared ratio to stay in exact rational arithmetic for every r.
        ratio_squared = Fraction(torsion * torsion, order)
        ratios.append(ratio_squared)
        assert torsion == 2**r
        assert ratio_squared == Fraction(2 ** (2 * r), 2 ** (r * (r + 1) // 2))
    assert ratios[-1] < ratios[0]
    assert ratios[-1] < Fraction(1, 2**20)


def check_cyclic_one_dimensional_kernel_bound() -> None:
    # For cyclic groups Z/m, every fixed multiplication kernel has size gcd(n,m),
    # hence is bounded by n, matching the one-dimensional counting mechanism.
    for n in range(1, 13):
        for m in range(1, 80):
            # Z/m is represented directly by its one cyclic factor, with prime
            # factorization supplied below.
            rem = m
            factors: list[Factor] = []
            p = 2
            while p * p <= rem:
                if rem % p == 0:
                    e = 0
                    while rem % p == 0:
                        rem //= p
                        e += 1
                    factors.append((p, e, 1))
                p += 1
            if rem > 1:
                factors.append((rem, 1, 1))
            assert kernel_size(n, factors) == __import__("math").gcd(n, m)
            assert kernel_size(n, factors) <= n


if __name__ == "__main__":
    check_decomposition_counts()
    check_bounded_width_obstruction()
    check_cyclic_one_dimensional_kernel_bound()
    print("decomposition_and_kernel_checks: OK")
    print("bounded_width_obstruction_checks: OK")
    print("cyclic_kernel_bound_checks: OK")
    print("ALL_CHECKS_OK")
