#!/usr/bin/env python3
"""Exact finite checks for Ballot--Kimberling--Moses, Table 5, (x-1)^2."""
from math import comb

MAX_DEPTH = 100
UNIVERSAL_DEPTH = 14
SOURCE_COUNTS = [1, 1, 2, 4, 7, 11, 16, 22, 29, 37, 46]


def normalize(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add_one(p):
    return normalize((p[0] + 1,) + p[1:])


def multiply_x(p):
    return normalize((0,) + p)


def remainder(p):
    """Integer long division by the monic polynomial X^2 - 2X + 1."""
    p = list(p)
    for j in range(len(p) - 1, 1, -1):
        c = p[j]
        p[j] = 0
        p[j - 1] += 2 * c
        p[j - 2] -= c
    p += [0, 0]
    return p[0], p[1]


def polynomial_fidelity():
    """Reproduce Theorem 1.1(1),(4) on complete small generations of T*."""
    frontier = {(0,)}
    visited = set(frontier)
    counts = [1]
    checked = 0
    for n in range(1, UNIVERSAL_DEPTH + 1):
        candidates = {f(p) for p in frontier for f in (add_one, multiply_x)}
        frontier = candidates - visited
        visited |= frontier
        counts.append(len(frontier))
        assert len(frontier) == 2 ** (n - 1), (n, len(frontier))
        for p in frontier:
            a = sum(p)
            b = sum(j * c for j, c in enumerate(p))
            assert n == a + len(p) - 1, (n, p)
            u, v = remainder(p)
            assert (u + v, v) == (a, b), (p, u, v, a, b)
            checked += 1
    print(f'FIDELITY T*: complete generations 0..{UNIVERSAL_DEPTH}; counts={counts}')
    print(f'FIDELITY T*: degree/weight and remainder-coordinate checks={checked}')


def check_tree():
    native_frontier = {(0, 0)}
    pair_frontier = {(0, 0)}
    native_seen = set(native_frontier)
    pair_seen = set(pair_frontier)
    counts = []
    distance_checks = 0
    for n in range(MAX_DEPTH + 1):
        translated = {(u + v, v) for u, v in native_frontier}
        assert translated == pair_frontier, ('two BFS models disagree', n)
        assert len(native_frontier) == 1 + comb(n, 2), (n, len(native_frontier))
        counts.append(len(native_frontier))
        for a, b in pair_frontier:
            if a == 0:
                assert (n, b) == (0, 0), (n, a, b)
                d = 0
            else:
                assert a >= 1 and b >= 0, (n, a, b)
                d = a + (b + a - 1) // a
            assert d == n, (n, a, b, d)
            distance_checks += 1
        if n == MAX_DEPTH:
            break
        native_candidates = {
            z for u, v in native_frontier
            for z in (remainder((u + 1, v)), remainder((0, u, v)))
        }
        pair_candidates = {
            z for a, b in pair_frontier for z in ((a + 1, b), (a, b + a))
        }
        native_frontier = native_candidates - native_seen
        pair_frontier = pair_candidates - pair_seen
        native_seen |= native_frontier
        pair_seen |= pair_frontier
    assert counts[:len(SOURCE_COUNTS)] == SOURCE_COUNTS
    for n in range(3, MAX_DEPTH + 1):
        assert counts[n] == 3 * counts[n-1] - 3 * counts[n-2] + counts[n-3], n
    assert len(native_seen) == len(pair_seen) == 166751
    print(f'FIDELITY Table 5, q=X^2-2X+1: G_0..G_10={counts[:11]}')
    print(f'CHECK quotient-remainder and dual-number BFS: complete depths 0..{MAX_DEPTH}')
    print(f'CHECK all observed first-generation distances: points={distance_checks}; mismatches=0')
    print(f'CHECK G_n=1+binomial(n,2): cases={MAX_DEPTH+1}; mismatches=0')
    print(f'CHECK recurrence: n=3..{MAX_DEPTH}; mismatches=0')
    print(f'COUNTS G_0..G_{MAX_DEPTH}={counts}')


def main():
    polynomial_fidelity()
    check_tree()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
