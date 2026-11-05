#!/usr/bin/env python3
"""Exact small-instance checks for Bowden--Tymoczko Question 6.3.

A divisor d of m denotes the ideal (d) in Z/mZ.  Thus d=m is the zero
ideal and d=1 is the unit ideal.  For a fixed labeling, the known CRT
rank formula is max_i c(A_i), where A_i contains edges whose label is
zero in Z/p_i^e_i Z.  The script independently enumerates every allowed
labeling for the requested small graphs and compares its minimum with the
partition formula.  It also enumerates the spline group for the source's
C3 rank-one example and computes its generator number from p-torsion.
"""

from __future__ import annotations

from collections import Counter
from itertools import product


def prime_powers(m: int) -> list[tuple[int, int, int]]:
    """Return (p, exponent, p**exponent) for the distinct prime factors of m."""
    out: list[tuple[int, int, int]] = []
    x = m
    p = 2
    while p * p <= x:
        if x % p == 0:
            q = 1
            e = 0
            while x % p == 0:
                x //= p
                q *= p
                e += 1
            out.append((p, e, q))
        p += 1
    if x > 1:
        out.append((x, 1, x))
    return out


def divisors(m: int) -> list[int]:
    return [d for d in range(1, m + 1) if m % d == 0]


def components(n: int, edges: list[tuple[int, int]], active: set[int]) -> int:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for j in active:
        u, v = edges[j]
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[rv] = ru
    return len({find(v) for v in range(n)})


def fixed_label_rank(n: int, edges: list[tuple[int, int]], labels: tuple[int, ...], m: int) -> int:
    """The CRT component formula for the additive generator number."""
    ranks = []
    for _, _, q in prime_powers(m):
        zero_edges = {j for j, d in enumerate(labels) if d % q == 0}
        ranks.append(components(n, edges, zero_edges))
    return max(ranks)


def canonical_partitions(n: int) -> list[tuple[frozenset[int], ...]]:
    """A canonical recursive partition generator (no duplicate block orderings)."""
    if n == 0:
        return [()]
    out: list[tuple[frozenset[int], ...]] = []
    for p in canonical_partitions(n - 1):
        out.append((frozenset({n - 1}),) + p)
        for i in range(len(p)):
            q = list(p)
            q[i] = q[i] | {n - 1}
            out.append(tuple(q))
    # The recursion above keeps the first-created block distinguished; sort and dedupe.
    seen = set()
    result = []
    for p in out:
        key = tuple(sorted((tuple(sorted(b)) for b in p)))
        if key not in seen:
            seen.add(key)
            result.append(tuple(frozenset(b) for b in key))
    return result


def crossing_edges(n: int, edges: list[tuple[int, int]], partition: tuple[frozenset[int], ...]) -> int:
    block = {v: i for i, b in enumerate(partition) for v in b}
    return sum(block[u] != block[v] for u, v in edges)


def partition_formula(n: int, edges: list[tuple[int, int]], m: int) -> int:
    t = len(prime_powers(m))
    values = []
    for partition in canonical_partitions(n):
        e = crossing_edges(n, edges, partition)
        values.append(len(partition) - ((t - 1) * e) // t)
    return max(values)


def all_splines(n: int, edges: list[tuple[int, int]], labels: tuple[int, ...], m: int) -> list[tuple[int, ...]]:
    """Enumerate the finite spline group for one small labeling exactly."""
    values = []
    for f in product(range(m), repeat=n):
        if all((f[u] - f[v]) % m % labels[j] == 0 for j, (u, v) in enumerate(edges)):
            values.append(f)
    return values


def direct_generator_number(group: list[tuple[int, ...]], m: int) -> int:
    """For a finite abelian group, max_p dim_Fp(G[p]) is d(G)."""
    if not group:
        return 0
    n = len(group[0])
    answers = []
    for p, _, _ in prime_powers(m):
        torsion = sum(all((p * x) % m == 0 for x in f) for f in group)
        # torsion is a power of p; use exact integer arithmetic only.
        q, exponent = 1, 0
        while q < torsion:
            q *= p
            exponent += 1
        assert q == torsion, (p, torsion)
        answers.append(exponent)
    return max(answers)


def graph_family() -> dict[str, tuple[int, list[tuple[int, int]]]]:
    out: dict[str, tuple[int, list[tuple[int, int]]]] = {}
    out["triangle"] = (3, [(0, 1), (1, 2), (0, 2)])
    out["K4"] = (4, [(i, j) for i in range(4) for j in range(i + 1, 4)])
    for n in range(2, 6):
        out[f"path{n}"] = (n, [(i, i + 1) for i in range(n - 1)])
    for n in range(3, 6):
        out[f"cycle{n}"] = (n, [(i, (i + 1) % n) for i in range(n)])
    return out


def enumerate_minima(n: int, edges: list[tuple[int, int]], m: int, allowed: list[int]):
    minimum = None
    witness = None
    distribution: Counter[int] = Counter()
    for labels in product(allowed, repeat=len(edges)):
        rank = fixed_label_rank(n, edges, labels, m)
        distribution[rank] += 1
        if minimum is None or rank < minimum:
            minimum, witness = rank, labels
    assert minimum is not None
    return minimum, witness, distribution


def main() -> None:
    ms = (4, 6, 12, 30)
    readings: list[str] = []

    # Source fidelity: Bowden--Tymoczko's C3 example with m = 2*3*5
    # and labels pq, qr, rp = 6, 15, 10 has rank one.
    c3_edges = [(0, 1), (1, 2), (0, 2)]
    example_labels = (6, 15, 10)
    example_group = all_splines(3, c3_edges, example_labels, 30)
    assert direct_generator_number(example_group, 30) == 1
    readings.append(f"fidelity C3,m=30,labels={example_labels}: rank=1, |R_G,alpha|={len(example_group)}")

    for name, (n, edges) in graph_family().items():
        for m in ms:
            ds = divisors(m)
            nonzero = [d for d in ds if d < m]
            proper = [d for d in nonzero if d > 1]
            unrestricted_min, unrestricted_witness, _ = enumerate_minima(n, edges, m, ds)
            nonzero_min, nonzero_witness, nonzero_dist = enumerate_minima(n, edges, m, nonzero)
            proper_min, proper_witness, proper_dist = enumerate_minima(n, edges, m, proper)
            # Every listed graph is connected, so the unrestricted all-ideal minimum is 1.
            assert unrestricted_min == 1, (name, m, unrestricted_min, unrestricted_witness)
            expected = partition_formula(n, edges, m)
            assert nonzero_min == expected, (name, m, nonzero_min, expected, nonzero_witness)
            assert proper_min == expected, (name, m, proper_min, expected, proper_witness)
            readings.append(
                f"{name:7s} m={m:2d}: unrestricted={unrestricted_min}; "
                f"nonzero={nonzero_min}; proper_nonzero={proper_min}; "
                f"formula={expected}; distribution_nonzero={dict(sorted(nonzero_dist.items()))}; "
                f"distribution_proper={dict(sorted(proper_dist.items()))}"
            )

    # The all-unit labeling is the unconstrained module (rank n), another fidelity check.
    for name, (n, edges) in graph_family().items():
        for m in ms:
            assert fixed_label_rank(n, edges, tuple(1 for _ in edges), m) == n

    print("Exact integer checks for Bowden--Tymoczko, arXiv:1501.02027v1")
    print("Convention tested: nonzero ideals d<m; proper nonzero ideals 1<d<m; unrestricted includes d=m.")
    for line in readings:
        print(line)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
