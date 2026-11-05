#!/usr/bin/env python3
"""Exact checks for Hamzeh--Ashrafi Conjecture 2.5 (IEJA 26, 2019).

No floating point and no third-party dependencies. A compressed cycle certificate
uses vertices (element order, fiber index); any enumeration of each actual order
fiber transports it to the literal graph on group elements.
"""
from collections import Counter, defaultdict
from itertools import permutations
from math import factorial, gcd, isqrt, lcm, prod

WORD = (1, 3, 6, 2, 4, 8, 4, 2, 6, 3, 9, 3, 15, 5, 10, 5,
        15, 3, 21, 7, 21, 3, 12)
EXPECTED_A10 = {1: 1, 2: 5355, 3: 31040, 4: 94500, 5: 78624,
                6: 201600, 7: 86400, 8: 226800, 9: 403200,
                10: 90720, 12: 302400, 15: 120960, 21: 172800}


def partitions(n, maximum=None):
    """Every integer partition, exactly once, including all fixed points."""
    if n == 0:
        yield ()
        return
    for first in range(min(n, n if maximum is None else maximum), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def fibers(n):
    counts = Counter()
    types = defaultdict(list)
    for part in partitions(n):
        if (n - len(part)) % 2:
            continue
        multiplicities = Counter(part)
        denominator = prod(j ** a * factorial(a)
                           for j, a in multiplicities.items())
        numerator = factorial(n)
        assert denominator > 0 and numerator % denominator == 0
        size = numerator // denominator
        m = lcm(*part)
        counts[m] += size
        types[m].append((part, size))
    assert sum(counts.values()) == factorial(n) // 2
    return dict(sorted(counts.items())), types


def comparable(a, b):
    assert a > 0 and b > 0
    return a % b == 0 or b % a == 0


def components(labels):
    remaining = set(labels)
    result = []
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        component, todo = {root}, [root]
        while todo:
            x = todo.pop()
            adjacent = {y for y in remaining if comparable(x, y)}
            remaining -= adjacent
            component |= adjacent
            todo.extend(adjacent)
        result.append(sorted(component))
    return result


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def has_label_hamilton_cycle(labels):
    """Exact subset dynamic programming; each label may be used only once."""
    labels = sorted(labels)
    n = len(labels)
    reachable = {1: {0}}
    for mask in range(1, 1 << n, 2):
        for last in tuple(reachable.get(mask, ())):
            for nxt in range(1, n):
                if not mask & (1 << nxt) and comparable(labels[last], labels[nxt]):
                    reachable.setdefault(mask | (1 << nxt), set()).add(nxt)
    return any(comparable(labels[last], labels[0])
               for last in reachable.get((1 << n) - 1, ()))


def certificate(counts, word):
    """Check block coverage, every expanded edge, and cyclic closure exactly."""
    multiplicity = Counter(word)
    assert set(multiplicity) == set(counts)
    occurrence = Counter()
    intervals = defaultdict(list)
    blocks = []
    for m in word:
        occurrence[m] += 1
        c, total, j = multiplicity[m], counts[m], occurrence[m]
        assert 1 <= c <= total
        lo, hi = ((1, total - c + 1) if j == 1
                  else (total - c + j, total - c + j))
        assert 1 <= lo <= hi <= total
        intervals[m].append((lo, hi))
        blocks.append((m, lo, hi))
    for m, pieces in intervals.items():
        cursor = 1
        for lo, hi in sorted(pieces):
            assert lo == cursor
            cursor = hi + 1
        assert cursor == counts[m] + 1
    previous = first = None
    vertex_count = 0
    for m, lo, hi in blocks:
        for i in range(lo, hi + 1):
            vertex = (m, i)
            if previous is None:
                first = vertex
            else:
                assert vertex != previous and comparable(previous[0], m)
            previous = vertex
            vertex_count += 1
    assert previous != first and comparable(previous[0], first[0])
    assert vertex_count == sum(counts.values())
    return multiplicity, vertex_count


def capacity_tree(counts):
    """The candidate's nonidentity order-label spanning tree, checked on data."""
    labels = set(counts) - {1}
    primes = sorted(m for m in labels if prime(m))
    edges = {(3, 6), (2, 6)}
    used = {2, 3, 6}
    for p in primes:
        if p >= 5:
            assert 3 * p in labels
            edges |= {(3, 3 * p), (p, 3 * p)}
            used |= {p, 3 * p}
    for m in sorted(labels - used):
        divisor = next(p for p in primes if m % p == 0)
        edges.add((divisor, m))
    assert len(edges) == len(labels) - 1
    degree = Counter()
    adjacent = defaultdict(set)
    for a, b in edges:
        assert a in labels and b in labels and a != b and comparable(a, b)
        degree.update((a, b))
        adjacent[a].add(b)
        adjacent[b].add(a)
    visited, todo = set(), [min(labels)]
    while todo:
        a = todo.pop()
        if a not in visited:
            visited.add(a)
            todo.extend(adjacent[a] - visited)
    assert visited == labels
    assert all(degree[m] <= counts[m] for m in labels)
    assert all(degree[m] <= 2 for m in labels if not prime(m))
    return dict(sorted(degree.items()))


def permutation_type(p):
    seen, lengths = set(), []
    for i in range(len(p)):
        if i not in seen:
            j, length = i, 0
            while j not in seen:
                seen.add(j)
                length += 1
                j = p[j]
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def fidelity():
    # Source Theorem 2.1 on cyclic groups: every degree is even iff |G| is odd.
    for n in range(3, 33):
        orders = [n // gcd(n, a) for a in range(n)]
        degrees = [sum(a != b and comparable(orders[a], orders[b])
                       for b in range(n)) for a in range(n)]
        assert all(d % 2 == 0 for d in degrees) == (n % 2 == 1)
    print('FIDELITY source Theorem 2.1: cyclic groups C_n, 3 <= n <= 32, PASS')
    # Source Theorem 2.3, instantiated on C6: a quotient cycle lifts.
    orders = [6 // gcd(6, a) for a in range(6)]
    quotient_cycle = (1, 2, 6, 3)
    cycle = [a for m in quotient_cycle for a in range(6) if orders[a] == m]
    assert sorted(cycle) == list(range(6))
    assert all(a != b and comparable(orders[a], orders[b])
               for a, b in zip(cycle, cycle[1:] + cycle[:1]))
    print('FIDELITY source Theorem 2.3: C6 quotient (1,2,6,3), lifted cycle', cycle)
    # Directly enumerate A5 to independently check cycle-type counting and 2.2.
    brute = Counter()
    for p in permutations(range(5)):
        part = permutation_type(p)
        if (5 - len(part)) % 2 == 0:
            brute[lcm(*part)] += 1
    formula, _ = fibers(5)
    assert dict(sorted(brute.items())) == formula == {1: 1, 2: 15, 3: 20, 5: 24}
    assert components(set(formula) - {1}) == [[2], [3], [5]]
    assert all(p * q not in formula for p, q in ((2, 3), (2, 5), (3, 5)))
    print('FIDELITY source Theorem 2.2: A5 has disconnected prime graph and',
          'three components after removing identity; order counts', formula)


def main():
    fidelity()
    counts, types = fibers(10)
    assert counts == EXPECTED_A10
    print('A10 order fibers (order: size; full even cycle partitions and sizes):')
    for m in counts:
        print(m, counts[m], types[m])
    assert len(counts) == 13 and sum(counts.values()) == 1814400
    assert list(counts) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 21]
    print('PRIOR FIDELITY arXiv:1707.07251v1 Theorem 2.31: A10 spectrum PASS')
    multiplicity, vertex_count = certificate(counts, WORD)
    print('A10 cyclic order word:', WORD)
    print('A10 word multiplicities:', dict(sorted(multiplicity.items())))
    print('A10 capacities (order, available, required):',
          [(m, counts[m], multiplicity[m]) for m in counts])
    print('A10 expanded certificate: vertices = edges =', vertex_count,
          '; all fibers covered once; every edge and closing edge checked')
    degree = capacity_tree(counts)
    print('A10 spanning-tree degrees:', degree)
    print('A10 nonidentity comparability components:', components(set(counts) - {1}))
    assert not has_label_hamilton_cycle(counts)
    print('A10 quotient is non-Hamiltonian: exact subset DP PASS')
    # Independent literal permutations for the error in the later isolation step.
    p = (1, 2, 3, 4, 5, 6, 0, 7, 8, 9)
    q = (1, 2, 3, 4, 5, 6, 0, 8, 9, 7)
    assert lcm(*permutation_type(p)) == 7
    assert lcm(*permutation_type(q)) == 21
    assert (10 - len(permutation_type(q))) % 2 == 0
    assert comparable(7, 21)
    print('LATER-THEOREM CHECK: (1 2 3 4 5 6 7)(8 9 10) is even, order 21, adjacent to order 7')
    print('BOUNDED TRIAGE: n, |orders|, components after identity deletion, capacity criterion')
    positive = []
    for n in range(5, 35):
        c, _ = fibers(n)
        forbidden = any(prime(t) for t in (n - 2, n - 1, n))
        connected = len(components(set(c) - {1})) == 1
        assert connected == (not forbidden)
        if not forbidden:
            d = capacity_tree(c)
            assert all(c[p] >= 2 ** (n - 1) for p in c if prime(p))
            assert all(c[m] >= d[m] for m in d)
            positive.append(n)
        print(n, len(c), len(components(set(c) - {1})),
              'PASS' if not forbidden else 'identity-cut obstruction')
    assert [6 * k + 4 for k in range(1, 6)] == [10, 16, 22, 28, 34]
    assert all(6 * k + 4 in positive for k in range(1, 6))
    print('Family prefix checked: A_(6k+4), 1 <= k <= 5; no unbounded computational claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
