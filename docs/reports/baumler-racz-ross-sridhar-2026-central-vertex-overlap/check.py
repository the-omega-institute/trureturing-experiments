#!/usr/bin/env python3
"""Exact finite tests; no finite test certifies the asymptotic conjecture.

Original standard-library code. Source fidelity: arXiv:2606.02472v1,
Definition 1.4 and the Pólya representation in the proof of Lemma 2.6.
All evidence arithmetic uses integers or fractions.Fraction.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def urn_pmf(n, k):
    law = {1: Q(1)}
    for time in range(k, n):
        nxt = Counter()
        for red, mass in law.items():
            nxt[red + 1] += mass * Q(red, time + 1)
            nxt[red] += mass * Q(time + 1 - red, time + 1)
        law = nxt
    return {s: p for s, p in law.items() if p}


def deleted_components(parents):
    N = len(parents)
    adj = [set() for _ in range(N)]
    for v in range(1, N):
        adj[v].add(parents[v]); adj[parents[v]].add(v)
    result = []
    for k in range(N):
        unseen, largest = set(range(N)), 0
        while unseen:
            seed = unseen.pop()
            seen, todo = {seed}, [seed]
            while todo:
                v = todo.pop()
                for w in adj[v] if v != k else ():
                    if w != k and w not in seen:
                        seen.add(w); unseen.discard(w); todo.append(w)
            largest = max(largest, len(seen))
        result.append(largest)
    return result


def overlap_bounds(Psi, K, r):
    """Extremal overlap across every possible cutoff tie selection."""
    if r == 0:
        return 0, 0
    cutoff = sorted(Psi, reverse=True)[r - 1]
    higher = [v for v, x in enumerate(Psi) if x > cutoff]
    tied = [v for v, x in enumerate(Psi) if x == cutoff]
    slots = r - len(higher)
    H = sum(v < K for v in higher)
    E = sum(v < K for v in tied)
    return H + max(0, slots - (len(tied) - E)), H + min(slots, E)


def main():
    total = covariance_checks = rank_checks = fidelity_checks = 0
    for n in range(9):
        N, count = n + 1, factorial(n)
        pmfs = [Counter() for _ in range(N)]
        hits = [[0] * N for _ in range(N + 1)]
        pairs = [[[0] * N for _ in range(N)] for _ in range(N + 1)]
        discrepancies = trees = 0
        for choices in product(*(range(j) for j in range(1, N))):
            parents = (0,) + choices
            F = [1] * N
            for v in range(n, 0, -1):
                F[parents[v]] += F[v]
            psi = [max([1, N - F[v]] +
                       [F[w] for w in range(1, N) if parents[w] == v])
                   for v in range(N)]
            if n <= 6:
                assert psi == deleted_components(parents)
            Psi = [N - x for x in psi]
            D = sum(x != y for x, y in zip(F, Psi))
            assert D <= sum(2 * x >= N for x in F)
            discrepancies += D
            for k in range(N):
                pmfs[k][F[k]] += 1
                if 2 * F[k] < N:
                    assert Psi[k] == F[k]
            for m in range(1, N + 1):
                active = [k for k in range(N) if F[k] >= m]
                for k in active:
                    hits[m][k] += 1
                for i, k in enumerate(active):
                    for ell in active[i + 1:]:
                        pairs[m][k][ell] += 1
            if n <= 6:
                for K in range(1, N + 1):
                    for r in range(N + 1):
                        lo, hi = overlap_bounds(Psi, K, r)
                        for m in range(1, N + 1):
                            q = sum(x >= m for x in Psi)
                            early = sum(Psi[v] >= m for v in range(K))
                            if q <= r:
                                assert early <= lo
                            if r <= q:
                                assert hi <= early
                            rank_checks += 1
            trees += 1
        assert trees == count
        total += trees
        assert Q(discrepancies, count) <= 2
        for k in range(1, N):
            assert {s: Q(v, count) for s, v in pmfs[k].items()} == urn_pmf(n, k)
            fidelity_checks += 1
        for m in range(1, N + 1):
            ps = [Q(hits[m][k], count) for k in range(N)]
            assert ps == [Q(choose(n - m + 1, k), choose(n, k)) for k in range(N)]
            assert all(ps[k] == Q(choose(n - k, m - 1), choose(n, m - 1))
                       for k in range(N))
            assert sum(ps) == Q(N, m)
            for k in range(1, N):
                for ell in range(k + 1, N):
                    cov = Q(pairs[m][k][ell], count) - ps[k] * ps[ell]
                    assert cov <= ps[ell] / (k + 1)
                    covariance_checks += 1
            for L in range(N):
                second = sum(ps[:L + 1]) + 2 * sum(
                    Q(pairs[m][k][ell], count)
                    for k in range(L + 1) for ell in range(k + 1, L + 1))
                variance = second - sum(ps[:L + 1]) ** 2
                bound = sum(ps[1:L + 1]) + 2 * sum(
                    ps[ell] / (k + 1)
                    for k in range(1, L + 1) for ell in range(k + 1, L + 1))
                assert 0 <= variance <= bound
        print(f'n={n}: maps={count}; tails/means/covariances/variance OK; '
              f'E[D]={Q(discrepancies, count)}')
    for K in [10, 100, 1000]:
        n, m = K * K, (K + 1) // 2
        mean = Q(sum(choose(n - k, m - 1) for k in range(K)), K * choose(n, m - 1))
        partial = Q(0)
        for j in range(41):
            partial += Q((-1) ** j, 2 ** j * factorial(j))
            if j == 39:
                lower_exp = partial
        lower_f, upper_f = 2 * (1 - partial), 2 * (1 - lower_exp)
        error = max(abs(mean - lower_f), abs(mean - upper_f))
        tolerance = Q(1, 20) if K == 10 else Q(1, 100)
        assert error < tolerance
        print(f'K={K}, n={n}, m={m}: E[S_n]/K={Q(n+1,m*K)}; '
              f'|E[S_(K-1)]/K-f(2)| < {tolerance} (exact rational enclosure)')
    print(f'Fidelity urn laws={fidelity_checks}; parent maps={total}; '
          f'covariance inequalities={covariance_checks}; tie-robust rank checks={rank_checks}')
    print('Finite checks support the route; they do not prove the asymptotic claim.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
