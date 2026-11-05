#!/usr/bin/env python3
"""Exact two-component positive-scattering checks for Wang--Yuan §7.1.

For two independent component vectors in one mode, the source's activation
formula gives kappa = +infinity for disjoint supports, kappa = 0 when a
shared normalized coordinate functional is a vertex, and kappa = 1 otherwise.
All arithmetic here is Fraction-based; no floating point is used.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path

INF = None


def support(v):
    return {i for i, x in enumerate(v) if x != 0}


def mode_cost(a, b):
    """Exact kappa_j({1,2}) from the source's normalized-dual test."""
    sa, sb = support(a), support(b)
    common = sa & sb
    if not common:
        return INF
    ratios = [Fraction(a[i], a[i] + b[i]) for i in range(len(a)) if a[i] + b[i] != 0]
    lo, hi = min(ratios), max(ratios)
    if any(Fraction(a[i], a[i] + b[i]) in (lo, hi) for i in common):
        return 0
    return 1


def tau_pair(modes):
    """Two-component activation formula: inf if two modes are impossible."""
    costs = [mode_cost(a, b) for a, b in modes]
    if sum(c is INF for c in costs) > 1:
        return INF
    return sum(c for c in costs if c is not INF)


def kronecker(vectors):
    out = [1]
    for v in vectors:
        out = [x * y for x in out for y in v]
    return out


def grouped(modes, partition):
    return [
        (kronecker([modes[i][0] for i in block]),
         kronecker([modes[i][1] for i in block]))
        for block in partition
    ]


def cat(name, modes, partition, expected_original, expected_grouped):
    got_original = tau_pair(modes)
    got_grouped = tau_pair(grouped(modes, partition))
    assert got_original == expected_original, (name, "original", got_original)
    assert got_grouped == expected_grouped, (name, "grouped", got_grouped)
    assert len(partition) >= 3, (name, "fewer than three grouped modes")
    print(f"{name}: tau={got_original}, tau_grouped={got_grouped}, grouped_modes={len(partition)}")


# A is a finite positive-scattering mode: both components have an exclusive
# support coordinate, with one shared coordinate.
A = ([1, 1, 0], [0, 1, 1])
# C is zero-cost: all coordinates are shared and the common support contains
# the normalized-dual extrema.
C = ([1, 1], [1, 2])
# P and Q are each zero-cost, but their opposite one-sided supports make
# P tensor Q a positive-scattering mode after grouping.
P = ([1, 1], [1, 0])
Q = ([1, 0], [1, 1])

for k in (1, 2, 3):
    # Increase family: k opposite P/Q pairs, one A mode, and one C mode.
    inc_modes = []
    inc_partition = []
    for j in range(k):
        start = len(inc_modes)
        inc_modes.extend([P, Q])
        inc_partition.append(tuple(range(start, start + 2)))
    a_index = len(inc_modes)
    inc_modes.append(A)
    c_index = len(inc_modes)
    inc_modes.append(C)
    inc_partition.extend([(a_index,), (c_index,)])
    cat(f"increase k={k}", inc_modes, inc_partition, 1, k + 1)

    # Decrease family: k+2 A modes and one C mode.  Group k+1 A modes
    # together, leaving one A and C as separate grouped modes.
    dec_modes = [A for _ in range(k + 2)] + [C]
    dec_partition = [tuple(range(k + 1)), (k + 1,), (k + 2,)]
    cat(f"decrease k={k}", dec_modes, dec_partition, k + 2, 2)

print("ALL_CHECKS_OK")
