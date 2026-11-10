#!/usr/bin/env python3
"""Exact finite checks for Rognes, arXiv:2607.00757v1, Section 4.

The source-fidelity check is the sum approximation of Example 2.2 and its
symmetrization in Theorem 2.1. The settlement checks exercise the proposed
separating formula and the gap-splitting lemma. Finite tests do not prove the
unbounded first-order invariance assertion. No third-party code is used.
"""

from itertools import product
from math import lcm


def source_sum(x, y, z):
    return (0 < y and x <= z) or (y == 0 and x == z) or (y < 0 and z <= x)


def source_symmetric_sum(x, y, z):
    return source_sum(x, y, z) and source_sum(y, x, z)


def sigma_plus(x, y, z):
    return x < y or y < x or z == 2 * x


def modular_sum(p, x, y, z):
    """The literal finite disjunction of Definition 2.3."""
    return any(
        x % p == a and y % p == b and z % p == (a + b) % p
        for a in range(p) for b in range(p)
    )


def phi_precise(x, y, u, v):
    return 0 < x < y < u < v and x + x == v and y + y == u


def phi_coarse(p, sigma, x, y, u, v):
    return (0 < x < y < u < v
            and modular_sum(p, x, x, v) and sigma(x, x, v)
            and modular_sum(p, y, y, u) and sigma(y, y, u))


def threshold(p, r):
    h = 1
    for _ in range(r):
        h = 2 * h + p
    return h


def gap_equivalent(h, a, b):
    return a == b or (a >= h and b >= h)


def tuple_equivalent(p, h, xs, ys):
    assert len(xs) == len(ys)
    return (all(x % p == y % p for x, y in zip(xs, ys))
            and all((x < x1) == (y < y1) and (x == x1) == (y == y1)
                    for x, y in zip(xs, ys) for x1, y1 in zip(xs, ys))
            and all(gap_equivalent(h, xs[i + 1] - xs[i], ys[i + 1] - ys[i])
                    for i in range(len(xs) - 1)))


def split_signatures(p, h, length):
    """All one-point extensions, including equality with either endpoint.

    Distances are truncated only after reaching h. The residue records the
    offset from the left endpoint. Equal sets certify both back-and-forth
    directions for these two finite intervals.
    """
    return frozenset((k % p, min(k, h), min(length - k, h))
                     for k in range(length + 1))


def main():
    fidelity = 0
    expanded = 0
    ordinary = 0
    for x, y in product(range(-64, 65), repeat=2):
        assert source_sum(x, y, x + y)
        assert source_symmetric_sum(x, y, x + y)
        fidelity += 1
        assert sigma_plus(x, y, x + y)
        expanded += 1
        # Addition determines u and v, so this includes all possible exact
        # solutions with x,y in this range, with no bound imposed on u,v.
        assert not phi_precise(x, y, 2 * y, 2 * x)
        ordinary += 1
    print(f'SOURCE_FIDELITY: {fidelity} integer input pairs satisfy Example 2.2 sum and Theorem 2.1 sum')
    print(f'PRECISE_PHI: {ordinary} integer input pairs have no exact completion')
    for x, z in product(range(-64, 65), range(-128, 129)):
        assert sigma_plus(x, x, z) == (z == 2 * x)
        expanded += 1
    print(f'EXPANDED_APPROXIMATION: {expanded} soundness/diagonal checks')

    coarse_count = 0
    for p in range(2, 18):
        for k in range(1, 41):
            n = p * k
            v = (n, 2 * n, 5 * n, 6 * n)
            assert phi_coarse(p, source_sum, *v)
            assert phi_coarse(p, source_symmetric_sum, *v)
            assert not phi_coarse(p, sigma_plus, *v)
            coarse_count += 1
    # One realization for a finite collection of tests with different moduli.
    n = lcm(*range(2, 18)) * 101
    for p in range(2, 18):
        assert phi_coarse(p, source_symmetric_sum, n, 2*n, 5*n, 6*n)
    print(f'COARSE_SEPARATION: {coarse_count} family members; finite collection moduli 2..17, N={n}')

    interval_pairs = 0
    for p in range(2, 14):
        for j in range(5):
            h = threshold(p, j)
            next_h = threshold(p, j + 1)
            assert next_h == 2 * h + p
            bound = next_h + 4 * p + 5
            signatures = [split_signatures(p, h, g) for g in range(bound + 1)]
            for g in range(1, bound + 1):
                for other in range(g, bound + 1):
                    if g % p == other % p and gap_equivalent(next_h, g, other):
                        assert signatures[g] == signatures[other], (p, j, g, other)
                        interval_pairs += 1
    print(f'RANK_SPLITTING: {interval_pairs} exact finite interval pairs; p=2..13, j=0..4, lengths<=H_(j+1)+4p+5')

    transfers = 0
    for p in range(2, 14):
        for r in range(7):
            for m in [0, 1, 3, 10, 100]:
                h = threshold(p, r)
                n = p * (m + h + 1)
                anchors = tuple(range(m + 1))
                assert tuple_equivalent(p, h, anchors + (n, 2*n), anchors + (n, 6*n))
                assert tuple_equivalent(p, h, anchors + (2*n, 4*n), anchors + (2*n, 5*n))
                transfers += 2
    print(f'TRANSFER_CONFIGURATIONS: {transfers} exact tuple invariants for N=p(M+H_r+1)')
    print('SCOPE: finite checks only; arbitrary quantified approximations require the unbounded invariance proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
