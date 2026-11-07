#!/usr/bin/env python3
"""Exact diagnostics for inverse-power line equilibria; no infinite-theorem certificate.

Source mathematics: A. Ferudun, arXiv:2610.06944v1, Lemma 2.3.
Independent implementation using only the Python standard library.
"""
from fractions import Fraction as Q
from itertools import combinations, product


def sign(x):
    return (x > 0) - (x < 0)


def conductance(a, n, j, s):
    d = abs(a[j] - a[n])
    e = abs(a[j + 1] - a[n + 1])
    assert d > 0 and e > 0
    return (d ** (-s) - e ** (-s)) / (e - d) if e != d else s * d ** (-s - 1)


def cutoff(x, r):
    return min(Q(1), max(Q(0), Q(2) - abs(x) / r))


def energy(edges, v):
    return sum((c * (v[i] - v[j]) ** 2 for i, j, c in edges), Q(0))


def source_fidelity():
    count = 0
    # Both branches of the source's exact inverse-square conductance formula.
    for d, e in product(map(Q, range(1, 10)), repeat=2):
        c = (d ** -2 - e ** -2) / (e - d) if d != e else 2 * d ** -3
        assert c == (d + e) / (d * d * e * e)
        count += 1
    print(f'FIDELITY Lemma 2.3 inverse-square formula: {count} exact pairs')


def truncated_harmonic_identity():
    count = 0
    for length in range(2, 6):
        for raw in product(range(1, 4), repeat=length):
            gaps = tuple(map(Q, raw))
            a = [Q(-sum(raw), 2)]
            for g in gaps:
                a.append(a[-1] + g)
            for s in range(2, 7):
                for n in range(length):
                    # Reindexing needs matching windows: [0,length) and [1,length].
                    # A general truncated configuration is not an equilibrium.
                    left = sum((sign(a[j] - a[n]) * abs(a[j] - a[n]) ** (-s)
                                for j in range(length) if j != n), Q(0))
                    shifted = sum((sign(a[j + 1] - a[n + 1]) *
                                   abs(a[j + 1] - a[n + 1]) ** (-s)
                                   for j in range(length) if j != n), Q(0))
                    rhs = sum((conductance(a, n, j, s) * (gaps[j] - gaps[n])
                               for j in range(length) if j != n), Q(0))
                    assert left - shifted == rhs
                    transformed_residual = sum((conductance(a, n, j, s) * gaps[n] * gaps[j] *
                                                (1 / gaps[j] - 1 / gaps[n])
                                                for j in range(length) if j != n), Q(0))
                    assert transformed_residual == -rhs
                    count += 1
    a = tuple(map(Q, [0, 1, 3, 6, 7]))
    residual = sum((conductance(a, 1, j, 3) * ((a[j+1]-a[j]) - 2)
                    for j in range(4) if j != 1), Q(0))
    assert residual != 0
    print(f'TRUNCATION aligned force-difference and reciprocal residual: {count} exact cases')
    print(f'TRUNCATION non-equilibrium example a=(0,1,3,6,7), n=1, s=3: residual={residual}')


def arithmetic_progression_checks():
    count = 0
    for alpha, beta, s, n, radius in product(
            [Q(1, 3), Q(1), Q(5, 2)], [Q(-2), Q(0), Q(7, 3)],
            range(2, 9), range(-2, 3), range(1, 13)):
        x = alpha * n + beta
        net = sum((sign(alpha * (n + k) + beta - x) *
                   abs(alpha * (n + k) + beta - x) ** (-s)
                   for k in range(-radius, radius + 1) if k), Q(0))
        assert net == 0
        # Exact integral of t^(-s) from radius to infinity, a one-sided tail bound.
        tail_integral = alpha ** (-s) / ((s - 1) * Q(radius) ** (s - 1))
        assert tail_integral > 0
        count += 1
    print(f'EXPECTED SETTLEMENT arithmetic-progression symmetric-force checks: {count} exact cases')


def ground_state_identity():
    count = 0
    for n in range(1, 4):
        pairs = list(combinations(range(n), 2))
        for weights in product([Q(0), Q(1)], repeat=len(pairs)):
            edges = [(i, j, c) for (i, j), c in zip(pairs, weights)]
            for f in product(map(Q, [-1, 0, 1]), repeat=n):
                for w in product(map(Q, [0, 1, 2]), repeat=n):
                    lap = [sum((c * (w[i] - w[j]) if i == x else
                                c * (w[j] - w[i]) if j == x else Q(0)
                                for i, j, c in edges), Q(0)) for x in range(n)]
                    lhs = energy(edges, [f[i] * w[i] for i in range(n)])
                    potential = sum((f[i] ** 2 * w[i] * lap[i] for i in range(n)), Q(0))
                    gradient = sum((c * w[i] * w[j] * (f[i] - f[j]) ** 2
                                    for i, j, c in edges), Q(0))
                    assert lhs == potential + gradient
                    count += 1
    # One induced vertex coupled to a higher-value exterior vertex; local harmonicity.
    h = [Q(1), Q(1, 2), Q(2)]
    c01, c02, m = Q(2), Q(1), Q(3, 2)
    assert c01 * (h[0] - h[1]) + c02 * (h[0] - h[2]) == 0
    w = [m - h[0], m - h[1]]
    induced_lap = c01 * (w[0] - w[1])
    assert induced_lap == -c02 * (h[2] - h[0]) < 0
    print(f'EXPECTED SETTLEMENT induced ground-state identity: {count} exhaustive exact cases')
    print('EXPECTED SETTLEMENT local sublevel flux sign: exact check OK')


def restricted_near_energy():
    count = 0
    r = Q(3)  # All neighboring gaps in this finite family have ratio in [1/3,3].
    for raw in product(range(1, 4), repeat=4):
        gaps = tuple(map(Q, raw))
        a = [Q(-sum(raw), 2)]
        for g in gaps:
            a.append(a[-1] + g)
        for s, eps, radius in product(range(3, 7), map(Q, [1, 2, 3]), map(Q, range(1, 7))):
            active = [n for n in range(4) if gaps[n] >= eps]
            near = Q(0)
            for n, j in combinations(active, 2):
                d = a[j] - a[n]
                if d <= radius:
                    near += gaps[n] * gaps[j] * d ** (-s-1) * (cutoff(a[n], radius)-cutoff(a[j], radius)) ** 2
                # Integral comparison on one interval, p=s-1>1.
                p = s - 1
                integral = (d ** (1-p) - (d+gaps[j]) ** (1-p)) / (p-1)
                assert gaps[j] * d ** (-p) <= (1+r) ** p * integral
            bound = 7 * (1+r) ** (s-1) * eps ** (2-s) / ((s-2) * radius)
            assert near <= bound
            count += 1
    print(f'EXPECTED SETTLEMENT restricted near-energy and integral bounds: {count} exact cases')


def main():
    source_fidelity()
    truncated_harmonic_identity()
    arithmetic_progression_checks()
    ground_state_identity()
    restricted_near_energy()
    print('SCOPE: finite algebra and bounds only; no certificate of the infinite rigidity theorem')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
