#!/usr/bin/env python3
"""Exact finite checks for Sibirmovsky's OEIS A276436 zero-member conjecture.

Source: https://oeis.org/A276436, revision 35; coefficients attributed to
Yuriy Sibirmovsky (OEIS data, CC BY-SA 4.0). The algorithms are independent
integer implementations. Finite output is not a proof of the universal claim.
"""
from collections import defaultdict
from math import comb, isqrt

DEGREE = 100
PROFILE_LIMIT = 100_000
ZERO_SET = {0, 3, 5, 8}
# The complete DATA field of OEIS revision 35, indexed from zero.
OEIS = [
    0, 1, 1, 0, 2, 0, 1, 6, 0, 6, 6, 24, 15, 26, 48, 56, 240, 60,
    303, 504, 780, 1002, 1776, 3246, 3601, 7826, 7500, 18980, 26874,
    38130, 56196, 99360, 153636, 210084, 390348, 486420, 900428,
    1310118, 2064612, 3073008, 4825138, 7558008, 11428162, 18596886,
    26006031, 43625940, 65162736, 100027728, 152897710, 242895050,
]


def multiply(a, b):
    out = [0] * (DEGREE + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(DEGREE + 1 - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def power(a, k):
    out = [1] + [0] * DEGREE
    for _ in range(k):
        out = multiply(out, a)
    return out


def approximation(height):
    """Source f_(height+1), evaluated modulo X^(DEGREE+1)."""
    out = [0, 1] + [0] * (DEGREE - 1)
    for depth in reversed(range(height)):
        out = power(out, depth + 2)
        out[1] += 1
    return out


def profile_coefficients():
    """Independent level-profile DP, weighting child-slot choices by binomials.

    A state at level d stores (sum_{j<=d} j*t_j, t_d). Stopping a state
    sets every later t_j to zero. Advancing chooses the next internal slots.
    """
    out = [0] * (DEGREE + 1)
    out[1] = 1  # Single-leaf tree.
    states = {(1, 1): 1}  # Internal root: t_1=1.
    d = 1
    states_visited = 0
    while states:
        states_visited += len(states)
        next_states = defaultdict(int)
        for (used, t), weight in states.items():
            out[1 + used] += weight
            slots = (d + 1) * t
            max_next = min(slots, (DEGREE - 1 - used) // (d + 1))
            for t_next in range(1, max_next + 1):
                next_states[(used + (d + 1) * t_next, t_next)] += (
                    weight * comb(slots, t_next)
                )
        states = next_states
        d += 1
    return out, states_visited


def triangular(h):
    return h * (h + 1) // 2


def construct_profile(n):
    """Construct t_1,...,t_H, or report an exceptional leaf count."""
    if n == 0:
        return None
    h = (isqrt(8 * (n - 1) + 1) - 1) // 2
    r = n - 1 - triangular(h)
    assert 0 <= r <= h
    if h == 0:
        return []
    if r == 1:
        if h < 4:
            return None
        t = [1] * (h - 1)
        t[1] = 2  # t_2.
        t[h - 2] = 2  # t_(h-1), distinct from t_2.
        return t
    t = [1] * h
    if r >= 2:
        t[r - 1] = 2
    return t


def check_profile(n, t):
    if not t:
        assert n == 1
        return
    assert t[0] == 1
    assert all(isinstance(v, int) and v >= 0 for v in t)
    assert all(t[d] <= (d + 1) * t[d - 1] for d in range(1, len(t)))
    assert n == 1 + sum(d * v for d, v in enumerate(t, start=1))
    # Independent edge-count identity: edges = vertices - 1.
    internal = sum(t)
    child_slots = sum((d + 1) * v for d, v in enumerate(t, start=1))
    assert child_slots == n + internal - 1


def main():
    polys = [approximation(h) for h in range(17)]
    reference = polys[15]
    assert reference == polys[16]
    assert len(OEIS) == 50
    assert reference[:50] == OEIS
    # The original MSE post additionally displays a(50)=365185374 for f_10.
    assert reference[50] == 365185374
    assert polys[9][:51] == OEIS + [365185374]
    print('FIDELITY: all 50 OEIS coefficients and the 51-term MSE f_10 prefix match')
    for n in range(1, 17):
        through = min(DEGREE, n * (n - 1) // 2)
        assert polys[n - 1][:through + 1] == reference[:through + 1]
    print('FIDELITY: source stabilization bound holds for f_1 through f_16, degree <= 100')
    zeros = {n for n, a in enumerate(reference) if a == 0}
    assert zeros == ZERO_SET
    print('POLYNOMIAL: zeros on 0..100 =', sorted(zeros))
    counts, states_visited = profile_coefficients()
    assert counts == reference
    print('PROFILE_DP: exact binomial-weight counts match all 101 polynomial coefficients')
    print('PROFILE_DP: states visited =', states_visited)
    print('COEFFICIENT: a(100) =', reference[100])
    supported = 0
    rejected = []
    for n in range(PROFILE_LIMIT + 1):
        t = construct_profile(n)
        if t is None:
            assert n in ZERO_SET
            rejected.append(n)
        else:
            assert n not in ZERO_SET
            check_profile(n, t)
            supported += 1
    assert rejected == sorted(ZERO_SET)
    assert supported == PROFILE_LIMIT - 3
    print('CONSTRUCTION: valid profiles for all', supported, 'nonexceptional indices in 1..100000')
    print('CONSTRUCTION: rejected indices =', rejected)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
