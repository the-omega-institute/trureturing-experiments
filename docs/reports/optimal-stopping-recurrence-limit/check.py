#!/usr/bin/env python3
"""Exact sanity checks for arXiv:2205.08495v1, Conjecture 1.

Finite rational checks support the preregistered analytic counterfamily; they
are not a proof of any limit or of the conjecture's negation.
"""
from fractions import Fraction as Q
from itertools import permutations
from math import factorial, isqrt

import sympy as S


def floor_q(x):
    return x.numerator // x.denominator


def source_fidelity():
    """Independently enumerate secretary outcomes, then test source §3.1."""
    permutation_count = 0
    probability_count = 0
    for n in range(1, 9):
        wins = [0] * (n + 1)
        for perm in permutations(range(1, n + 1)):
            permutation_count += 1
            for k in range(n):
                best = max(perm[:k], default=0)
                selected = next((a for a in perm[k:] if a > best), None)
                wins[k] += selected == n
        probabilities = [Q(w, factorial(n)) for w in wins]
        assert probabilities[n] == 0
        for k in range(n):
            assert probabilities[k] == Q(1, n) + Q(k, k + 1) * probabilities[k + 1]
            probability_count += 1
        for k in range(1, n):
            assert probabilities[k] == Q(k, n) * sum(
                (Q(1, j) for j in range(k, n)), Q(0)
            )
    x = S.symbols('x', positive=True)
    f = -x * S.log(x)
    assert S.simplify(S.diff(f, x) - (f / x - 1)) == 0
    assert f.subs(x, 1) == 0
    print('SOURCE_FIDELITY: n=1..8; permutations=%d; recurrence checks=%d; '
          'f=-x*log(x) solves source ODE and f(1)=0' %
          (permutation_count, probability_count))


def make_row(n, alpha, theta, q, c):
    j = floor_q(alpha * n) + isqrt(n)
    active = 1 <= j < n
    H, G = [], []
    for k in range(n):
        spike = int(active and k == j)
        h = Q(n - k, n - k + 1) * (1 - q * spike)
        H.append(h)
        G.append(theta * (1 - h) + c * spike)
    F = [Q(0)] * n + [theta]
    for k in range(n - 1, -1, -1):
        F[k] = G[k] + H[k] * F[k + 1]
    return j, active, H, G, F


def family_checks():
    # Rational parameter tuples strictly inside all preregistered inequalities.
    instances = [
        (Q(1, 2), Q(1, 4), Q(1, 2), Q(1, 4), 500),
        (Q(1, 3), Q(2, 5), Q(2, 3), Q(1, 5), 160),
        (Q(3, 4), Q(1, 5), Q(3, 7), Q(1, 7), 160),
        (Q(1, 8), Q(4, 5), Q(1, 3), Q(1, 30), 160),
    ]
    total_checks = 0
    for alpha, theta, q, c, max_n in instances:
        assert all(0 < t < 1 for t in (alpha, theta, q))
        assert 0 < c < (1 - theta) * q
        checks, active_rows = 0, 0
        for n in range(1, max_n + 1):
            j, active, H, G, F = make_row(n, alpha, theta, q, c)
            assert F[n] == theta
            active_rows += active
            for k in range(n):
                assert 0 < H[k] < 1
                assert 0 < G[k] < 1 - H[k]
                assert 0 < G[k] / (1 - H[k]) < 1
                closed = theta + (c * Q(n - j + 1, n - k + 1)
                                  if active and k <= j else 0)
                assert F[k] == closed
                assert F[k] == G[k] + H[k] * F[k + 1]
                assert theta <= F[k] <= theta + c < 1
                assert n * G[k] - theta * n * (1 - H[k]) == (
                    n * c if active and k == j else 0)
                if not (active and k == j):
                    assert n * (1 - H[k]) == Q(n, n - k + 1)
                checks += 1
            if active:
                assert F[j] == theta + c
                assert all(F[k] < F[j] for k in range(n + 1) if k != j)
                assert floor_q(alpha * n) != j
            # Source Theorem 2(3iii), with the limiting solution f=theta:
            # g((k+1)/n) - theta*h((k+1)/n) = 0 exactly.
            residual_mass = sum((abs(n * G[k] - theta * n * (1 - H[k]))
                                 for k in range(1, n - 1)), Q(0)) / n
            assert residual_mass == (c if active and j <= n - 2 else 0)
        total_checks += checks
        print('FAMILY: alpha=%s theta=%s q=%s c=%s; n=1..%d; '
              'row-entry checks=%d; active rows=%d' %
              (alpha, theta, q, c, max_n, checks, active_rows))

    # Exact square subsequence: n=m^2, j=n/2+m. This also gives a
    # rigorously decaying rational error at x=alpha, rather than a decimal plot.
    alpha, theta, q, c = instances[0][:4]
    for m in (4, 8, 16, 32, 64):
        n = m * m
        k = n // 2
        j = k + m
        F_alpha = theta + c * Q(n - j + 1, n - k + 1)
        gap_to_left_limit = theta + c - F_alpha
        assert gap_to_left_limit == Q(m, 2 * (m * m + 2))
        assert 0 < gap_to_left_limit < Q(1, 2 * m)
        x_right = Q(3, 4)
        assert floor_q(n * x_right) >= j
        if m > 4:
            assert floor_q(n * x_right) > j
        print('EXACT_INTERIOR: n=%d; F_n(floor(n/2))=%s; '
              'gap to 1/2=%s < 1/(2m)' % (n, F_alpha, gap_to_left_limit))
    print('TOTAL_FAMILY_ROW_ENTRY_CHECKS=%d' % total_checks)


def symbolic_checks():
    x, theta, C, alpha, c = S.symbols('x theta C alpha c', real=True)
    y = theta + C / (1 - x)
    assert S.simplify(S.diff(y, x) - (y - theta) / (1 - x)) == 0
    assert S.simplify((1 - x) * (y - theta)) == C
    # Product-rule cancellation for an arbitrary ODE solution.
    Y = S.Function('Y')(x)
    invariant_derivative = S.diff((1 - x) * (Y - theta), x)
    assert S.simplify(invariant_derivative.subs(S.diff(Y, x),
                                               (Y - theta) / (1 - x))) == 0
    left_limit = theta + c * (1 - alpha) / (1 - x)
    assert S.simplify(left_limit.subs(x, alpha) - theta) == c
    assert S.limit((1 - x) * (y - theta), x, 1, dir='-') == C
    assert S.diff(theta, x) == 0
    print('SYMBOLIC: ODE invariant derivative=0; general solution verified; '
          'jump=c; endpoint constant solution verified')


def main():
    source_fidelity()
    family_checks()
    symbolic_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
