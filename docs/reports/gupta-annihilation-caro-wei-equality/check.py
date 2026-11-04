#!/usr/bin/env python3
"""Exact finite checks for Gupta's TxGraffiti equality classification.

Python standard library only. The exhaustive domain is every labelled simple
connected graph of order 3 through 6. Long path/cycle checks are separate
finite evidence; they do not prove the unbounded classification.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations


def graph(n, edges):
    adj = [0] * n
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return tuple(adj)


def path(n):
    return graph(n, ((i, i + 1) for i in range(n - 1)))


def cycle(n):
    assert n >= 3
    return graph(n, [(i, i + 1) for i in range(n - 1)] + [(n - 1, 0)])


def complete(n):
    return graph(n, combinations(range(n), 2))


def connected(adj):
    seen, todo = 0, 1
    while todo:
        bit = todo & -todo
        todo -= bit
        if seen & bit:
            continue
        seen |= bit
        todo |= adj[bit.bit_length() - 1] & ~seen
    return seen == (1 << len(adj)) - 1


def alpha_subsets(adj):
    # Independent-set recurrence, with each subset checked exactly once.
    independent = [True] * (1 << len(adj))
    best = 0
    for s in range(1, len(independent)):
        bit = s & -s
        rest = s ^ bit
        independent[s] = independent[rest] and not (adj[bit.bit_length() - 1] & rest)
        if independent[s]:
            best = max(best, s.bit_count())
    return best


def annihilation(degrees):
    m, total, answer = sum(degrees) // 2, 0, 0
    assert sum(degrees) % 2 == 0
    for d in sorted(degrees):
        if total + d > m:
            break
        total += d
        answer += 1
    return answer


def annihilation_subsets(degrees):
    # Independent implementation: maximize over every subset, not sorted prefixes.
    m = sum(degrees) // 2
    return max(s.bit_count() for s in range(1 << len(degrees))
               if sum(d for i, d in enumerate(degrees) if s >> i & 1) <= m)


def residue_list(degrees):
    ds = list(degrees)
    while True:
        ds.sort(reverse=True)
        if not ds or ds[0] == 0:
            return len(ds)
        d = ds.pop(0)
        assert d <= len(ds) and all(x >= 1 for x in ds[:d]), ds
        ds[:d] = [x - 1 for x in ds[:d]]


def residue_counts(degrees):
    # Independent implementation: act on multiplicities; delay inserting all
    # decremented entries until the entire original d-largest batch is chosen.
    counts = Counter(degrees)
    while counts:
        d = max(counts)
        if d == 0:
            return counts[0]
        counts[d] -= 1
        if counts[d] == 0:
            del counts[d]
        pending = Counter()
        for _ in range(d):
            k = max(counts)
            assert k > 0
            counts[k] -= 1
            if counts[k] == 0:
                del counts[k]
            pending[k - 1] += 1
        counts.update(pending)
    return 0


def invariants(adj, alpha=None):
    deg = [a.bit_count() for a in adj]
    D = max(deg, default=0)
    A = alpha_subsets(adj) if alpha is None else alpha
    a, R = annihilation(deg), residue_list(deg)
    assert R == residue_counts(deg)
    W = sum((Fraction(1, d + 1) for d in deg), Fraction())
    return D, A, a, R, W


def isomorphic(a, b):
    if len(a) != len(b):
        return False
    da, db = tuple(x.bit_count() for x in a), tuple(x.bit_count() for x in b)
    if sorted(da) != sorted(db):
        return False
    for p in permutations(range(len(a))):
        if any(da[i] != db[p[i]] for i in range(len(a))):
            continue
        if all(bool(a[i] >> j & 1) == bool(b[p[i]] >> p[j] & 1)
               for i in range(len(a)) for j in range(i)):
            return True
    return False


TARGETS = {'K4': complete(4), 'P3': path(3), 'P4': path(4), 'P6': path(6),
           'C3': cycle(3), 'C4': cycle(4), 'C5': cycle(5), 'C7': cycle(7)}


def matching_targets(adj):
    return [name for name, target in TARGETS.items() if isomorphic(adj, target)]


def binary_chain_alpha(n, closed):
    # Maximum-weight independent set on a binary chain, by two-state transfer;
    # closing the chain forbids simultaneous first and last selection.
    result = 0
    for first in [0, 1]:
        states = {first: first}
        for _ in range(1, n):
            nxt = {}
            for previous, score in states.items():
                for current in [0, 1]:
                    if not (previous and current):
                        nxt[current] = max(nxt.get(current, -1), score + current)
            states = nxt
        result = max(result, max(score for last, score in states.items()
                                 if not (closed and first and last)))
    return result


def main():
    print('Exact arithmetic: Python integers and fractions.Fraction; no floats.')
    k4 = invariants(complete(4))
    assert k4 == (3, 1, 2, 1, Fraction(1))
    print('FIDELITY Gupta Corollary 1, K4: (Delta,alpha,a,R,W) =', k4)
    k2 = invariants(complete(2))
    assert k2[:4] == (1, 1, 1, 1) and k2[0] * k2[1] < k2[2] + k2[3]
    print('FIDELITY Gupta Theorem 2 excluded K2: Delta*alpha=1 < a+R=2.')
    for name, adj in TARGETS.items():
        D, A, a, R, W = invariants(adj)
        assert D * A == a + R and matching_targets(adj) == [name]
        print('TARGET', name, '(Delta,alpha,a,R,W) =', (D, A, a, R, W))
    expected_counts = {3: 4, 4: 38, 5: 728, 6: 26704}
    total, equality_counts = 0, Counter()
    for n in range(3, 7):
        pairs = list(combinations(range(n), 2))
        tested = 0
        for mask in range(1 << len(pairs)):
            adj = graph(n, (e for i, e in enumerate(pairs) if mask >> i & 1))
            if not connected(adj):
                continue
            tested += 1
            D, A, a, R, W = invariants(adj)
            assert a == annihilation_subsets([x.bit_count() for x in adj])
            # Source Theorem 1, Theorem 2, and its three classical comparisons.
            assert a <= Fraction(D + 1, 2) * W
            assert W <= R <= A <= a
            assert D * A >= a + R
            # Classical Caro-Wei equality rigidity, separately checked.
            assert (W == A) == (adj == complete(n))
            names = matching_targets(adj)
            assert (D * A == a + R) == bool(names), (n, mask, D, A, a, R, names)
            equality_counts.update(names)
        assert tested == expected_counts[n]
        total += tested
        print('EXHAUSTIVE labelled connected n =', n, 'graphs =', tested)
    print('EXHAUSTIVE total =', total, 'equality occurrences =', dict(sorted(equality_counts.items())))
    paths, cycles = [], []
    for n in range(3, 301):
        for closed in [False, True]:
            adj = cycle(n) if closed else path(n)
            A = binary_chain_alpha(n, closed)
            D, _, a, R, W = invariants(adj, alpha=A)
            assert A == (n // 2 if closed else (n + 1) // 2) and a == A
            assert R == ((n + 2) // 3 if closed else n // 3 + 1)
            assert (D * A == a + R) == (n in ({3, 4, 5, 7} if closed else {3, 4, 6}))
            if D * A == a + R:
                (cycles if closed else paths).append(n)
    print('PATHS 3<=n<=300 equality orders =', paths)
    print('CYCLES 3<=n<=300 equality orders =', cycles)
    print('Residue: sorted-list and independent multiplicity implementations agree throughout.')
    print('Scope: finite evidence only; the preregistered settlement needs a general argument.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
