#!/usr/bin/env python3
"""Exact finite checks for Aprile--Ferrarini's strong total cover demand question.

All inequalities use integers after multiplication by a positive denominator.
Fraction is used for displayed LP and dual values. No floating-point solver or
third-party implementation is used. These finite checks are not a proof of the
unbounded rounding claim.
"""
from fractions import Fraction
from itertools import combinations, product
from random import Random


def graphs(n):
    potential = tuple(combinations(range(n), 2))
    for bits in product((0, 1), repeat=len(potential)):
        yield tuple(e for e, b in zip(potential, bits) if b)


def rows(n, edges):
    """Vertex variables precede edge variables; edges and vertices are disjoint."""
    result = [tuple([v] + [n + i for i, e in enumerate(edges) if v in e])
              for v in range(n)]
    result += [(u, v, n + i) for i, (u, v) in enumerate(edges)]
    assert all(len(r) == 3 and len(set(r)) == 3 for r in result[n:])
    assert all(sum(a in r for r in result[:n]) == (1 if a < n else 2)
               for a in range(n + len(edges)))
    return tuple(result)


def binary_round(n, rr, c, d, z, den):
    """z/den is boxed; c may be scaled by any common positive denominator."""
    assert all(0 <= a <= den for a in z)
    assert all(0 <= a for a in c)
    assert all(d[i] * den <= sum(z[a] for a in r) for i, r in enumerate(rr))
    high = {a for a, za in enumerate(z) if 3 * za >= den}
    # All feasible edge demands, including non-maximal ones, are covered by high.
    for i, r in enumerate(rr[n:], n):
        assert len(high.intersection(r)) >= d[i]
    selected = set(high)
    repair_cost = 0
    for v, r in enumerate(rr[:n]):
        deficit = max(0, d[v] - len(high.intersection(r)))
        low = sorted((a for a in r if a not in high), key=lambda a: (c[a], a))
        assert deficit * den <= sum(z[a] for a in low)
        assert deficit <= len(low)
        chosen = low[:deficit]
        assert sum(c[a] for a in chosen) * den <= sum(c[a] * z[a] for a in low)
        repair_cost += sum(c[a] for a in chosen)
        selected.update(chosen)
    high_fractional = sum(c[a] * z[a] for a in high)
    low_fractional = sum(c[a] * z[a] for a in range(len(z)) if a not in high)
    assert sum(c[a] for a in high) * den <= 3 * high_fractional
    assert repair_cost * den <= 2 * low_fractional
    x = tuple(int(a in selected) for a in range(len(z)))
    assert all(sum(x[a] for a in r) >= d[i] for i, r in enumerate(rr))
    assert sum(c[a] * x[a] for a in range(len(z))) * den <= 3 * sum(
        c[a] * z[a] for a in range(len(z)))
    return x


def integer_round(n, rr, c, d, z, den):
    assert all(0 <= a for a in z)
    assert all(d[i] * den <= sum(z[a] for a in r) for i, r in enumerate(rr))
    q = tuple(a // den for a in z)
    f = tuple(a % den for a in z)
    residual = tuple(max(0, d[i] - sum(q[a] for a in r)) for i, r in enumerate(rr))
    y = binary_round(n, rr, c, residual, f, den)
    x = tuple(a + b for a, b in zip(q, y))
    assert all(sum(x[a] for a in r) >= d[i] for i, r in enumerate(rr))
    assert sum(c[a] * x[a] for a in range(len(z))) * den <= 3 * sum(
        c[a] * z[a] for a in range(len(z)))
    return x


def fidelity():
    # Source Corollary (complete graphs): beta_T^*(K_N) = N-1, N >= 2.
    for n in range(2, 6):
        ee = tuple(combinations(range(n), 2))
        rr = rows(n, ee)
        k = n + len(ee)
        best = min(sum(x) for x in product((0, 1), repeat=k)
                   if all(sum(x[a] for a in r) >= 1 for r in rr))
        assert best == n - 1
        print(f'FIDELITY complete K_{n}: integer optimum={best}')
    # Source Proposition (WTC gap): K_{n+1}, zero-cost star at vertex 0.
    for n in range(2, 6):
        nv = n + 1
        ee = tuple(combinations(range(nv), 2))
        rr = rows(nv, ee)
        k = nv + len(ee)
        c = tuple([1] * nv + [int(0 not in e) for e in ee])
        zero = {a for a in range(k) if c[a] == 0}
        positive = tuple(a for a in range(k) if c[a])
        best = k
        for bits in product((0, 1), repeat=len(positive)):
            objective = sum(bits)
            if objective >= best:
                continue
            s = zero.union(a for a, b in zip(positive, bits) if b)
            if all(s.intersection(r) for r in rr):
                best = objective
        assert best == n - 1
        z = tuple([Fraction(0)] + [Fraction(1, 2)] * n +
                  [Fraction(int(0 in e)) for e in ee])
        assert all(sum(z[a] for a in r) >= 1 for r in rr)
        primal = sum(c[a] * z[a] for a in range(k))
        # A feasible exact dual independently certifies this LP optimum.
        dual = tuple([Fraction(0)] * nv +
                     [Fraction(int(0 not in e), n - 1) for e in ee])
        assert all(sum(dual[i] for i, r in enumerate(rr) if a in r) <= c[a]
                   for a in range(k))
        assert primal == sum(dual) == Fraction(n, 2)
        print(f'FIDELITY weighted K_{nv}: IP={best}, LP={primal}, ratio={Fraction(best)/primal}')


def exhaustive_binary():
    counts = {}
    for n in range(4):
        count = 0
        for ee in graphs(n):
            rr = rows(n, ee)
            k = len(rr)
            # All coordinates on the thirds grid and all 0/1 cost numerators.
            # Costs can be read as c/2, so nonintegral rational costs are included.
            for z in product(range(4), repeat=k):
                capacity = tuple(sum(z[a] for a in r) // 3 for r in rr)
                for vd in product(*(range(capacity[v] + 1) for v in range(n))):
                    # Edge demand does not enter rounding. Checking its greatest
                    # feasible value proves the same finite test for every smaller
                    # edge demand, without repeating an identical computation.
                    d = tuple(vd) + capacity[n:]
                    for c in product((0, 1), repeat=k):
                        binary_round(n, rr, c, d, z, 3)
                        count += 1
        counts[n] = count
        print(f'BINARY exhaustive n={n}: {count} cases; all vertex demands; edge demands by monotonicity')
    return sum(counts.values())


def exhaustive_integer():
    total = 0
    for n in range(3):
        count = 0
        for ee in graphs(n):
            rr = rows(n, ee)
            k = len(rr)
            # Every coordinate in {0, 1/3, ..., 8/3}; every feasible demand;
            # every cost numerator in {0,1}, again allowing the scale c/2.
            for z in product(range(9), repeat=k):
                capacity = tuple(sum(z[a] for a in r) // 3 for r in rr)
                for d in product(*(range(m + 1) for m in capacity)):
                    for c in product((0, 1), repeat=k):
                        integer_round(n, rr, c, d, z, 3)
                        count += 1
        total += count
        print(f'INTEGER exhaustive n={n}: {count} cases; every feasible demand vector')
    return total


def additional_checks():
    rng = Random(261023312)
    count = 0
    for _ in range(10000):
        n = rng.randrange(1, 9)
        ee = tuple(e for e in combinations(range(n), 2) if rng.randrange(3) == 0)
        rr = rows(n, ee)
        k = len(rr)
        z = tuple(rng.randrange(49) for _ in range(k))  # /12, includes 1/4 and 1/3
        c = tuple(rng.randrange(14) for _ in range(k))  # /7, exact rational costs
        capacity = tuple(sum(z[a] for a in r) // 12 for r in rr)
        d = tuple(rng.randrange(m + 1) for m in capacity)
        integer_round(n, rr, c, d, z, 12)
        boxed = tuple(a % 13 for a in z)
        cap = tuple(sum(boxed[a] for a in r) // 12 for r in rr)
        bd = tuple(rng.randrange(m + 1) for m in cap)
        binary_round(n, rr, c, bd, boxed, 12)
        count += 2
    print(f'ADDITIONAL seeded exact tests: {count}; n<=8; costs /7; coordinates /12')
    return count


if __name__ == '__main__':
    fidelity()
    total = exhaustive_binary() + exhaustive_integer() + additional_checks()
    print(f'ROUNDING_CHECKS={total}')
    print('ALL_CHECKS_OK')
