#!/usr/bin/env python3
"""Exact checks for the Fici--Shallit--Simpson asymptotic enumeration proposal.

Source: arXiv:2407.10564v2, section 7; CPM 2025, article 11, section 7.
No finite check establishes an asymptotic theorem.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product

SOURCE = [2, 4, 8, 16, 32, 58, 108, 190, 336, 560, 948, 1574,
          2568, 4116, 6596, 10444, 16320, 25488, 39216, 60690, 92204]
OEIS_EXTRA = [141060, 212104, 322508]


@lru_cache(None)
def palindromes(d):
    return tuple(s + s[:d // 2][::-1]
                 for t in product('01', repeat=(d + 1) // 2)
                 for s in [''.join(t)])


def words(n, m, a):
    return {(u * ((n + m - 1) // m))[:n]
            for p in palindromes(a) for s in palindromes(m - a)
            for u in [p + s]}


def literal_predicate(w):
    # Independent recognition: try every prefix length and every cut, including
    # the right-empty cut omitted by the presentation enumeration.
    n = len(w)
    for m in range(1, n + 1):
        u = w[:m]
        if any(w[i] != u[i % m] for i in range(n)):
            continue
        for a in range(m + 1):
            p, s = u[:a], u[a:]
            if p == p[::-1] and s == s[::-1]:
                return True
    return False


def matching(m, a):
    return frozenset((i, a - 1 - i) for i in range(a // 2)) | frozenset(
        (a + i, m - 1 - i) for i in range((m - a) // 2))


def constraints(n, m, a):
    return matching(m, a) | frozenset((i % m, i) for i in range(m, n))


def components(n, edges):
    parent = list(range(n))
    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for i, j in edges:
        i, j = find(i), find(j)
        parent[i] = j
    return len({find(i) for i in range(n)})


def T(m):
    return (3 * (m // 2) * 2 ** (m // 2) if m % 2 == 0
            else m * 2 ** ((m + 1) // 2))


def main():
    counts = []
    for n in range(1, 25):
        union = set()
        mass = 0
        for m in range(1, n + 1):
            for a in range(m):
                event = words(n, m, a)
                assert len(event) == 2 ** ((a + 1) // 2 + (m - a + 1) // 2)
                mass += len(event)
                union.update(event)
        assert mass == sum(T(m) for m in range(1, n + 1))
        if n % 2 == 0:
            k = n // 2
            assert mass == (10 * k - 12) * 2 ** k + 12
            lead = 10 * k * 2 ** k
        else:
            k = n // 2
            assert mass == (14 * k - 10) * 2 ** k + 12
            lead = 7 * (2 * k + 1) * 2 ** k
        assert len(union) <= mass
        counts.append(len(union))
        if n <= 14:
            brute = {''.join(t) for t in product('01', repeat=n)
                     if literal_predicate(''.join(t))}
            assert brute == union
        print(f'n={n} A={len(union)} presentation_mass={mass} A/leading={Fraction(len(union), lead)}')
    assert counts[:21] == SOURCE
    assert counts[21:] == OEIS_EXTRA
    print('FIDELITY_OK: all 21 arXiv source values and OEIS terms 22..24')
    print('LITERAL_CROSSCHECK_OK: all binary words of lengths 1..14')
    assert literal_predicate('010110')
    assert not literal_predicate('010110010110')
    assert not literal_predicate('001011')
    assert literal_predicate('0100110')
    print('SOURCE_EXAMPLES_OK')

    pair_checks = 0
    for n in range(1, 17):
        ps = [(m, a) for m in range(1, n + 1) for a in range(m)]
        for (m, a), (l, b) in combinations(ps, 2):
            M, L = matching(m, a), matching(l, b)
            e, S = len(M & L), len(M) + len(L)
            c = components(n, M | L)
            assert 4 * c <= 4 * n - 3 * S + 2 * e
            assert components(n, constraints(n, m, a) | constraints(n, l, b)) <= c
            if m == l:
                assert e == 0
            else:
                # Normalize the longer root to m.
                if m < l:
                    m, a, l, b = l, b, m, a
                delta = a - b
                assert e == 0 or delta in (0, l, l - m)
                if delta == 0:
                    assert e == a // 2
                if delta == l - m:
                    assert e == (l - b) // 2
            if n <= 8:
                left, right = words(n, m, a), words(n, l, b)
                assert len(left & right) == 2 ** components(
                    n, constraints(n, m, a) | constraints(n, l, b))
            pair_checks += 1
    print(f'MATCHING_AND_INTERSECTION_OK: {pair_checks} presentation pairs; explicit sets through n=8')

    cutoff_checks = 0
    for n in range(4, 25):
        for D in range(4):
            ps = [(m, a) for m in range(n - D, n + 1) for a in range(m)]
            near_mass = sum(T(m) for m in range(n - D, n + 1))
            pair_sum = sum(2 ** components(n, constraints(n, m, a) | constraints(n, l, b))
                           for (m, a), (l, b) in combinations(ps, 2))
            relaxed_sum = sum(2 ** components(n, matching(m, a) | matching(l, b))
                              for (m, a), (l, b) in combinations(ps, 2))
            # Integer relaxation of O_D(2^(n/2) + n^2 2^(n/4)).
            quarter = (n + 3) // 4
            cap = (D + 1) ** 2 * 2 ** (quarter + D + 2) * (
                n ** 2 + (32 + D) * 2 ** quarter)
            assert pair_sum <= relaxed_sum <= cap
            assert near_mass - pair_sum <= counts[n - 1]
            cutoff_checks += 1
    print(f'CUTOFF_BONFERRONI_OK: {cutoff_checks} cases, 4<=n<=24 and 0<=D<=3')
    for k in range(1, 101):
        assert sum(T(m) for m in range(1, 2 * k + 1)) == (10 * k - 12) * 2 ** k + 12
        assert sum(T(m) for m in range(1, 2 * k + 2)) == (14 * k - 10) * 2 ** k + 12
    print('PARITY_SUM_IDENTITIES_OK: 1<=k<=100')
    print('ASYMPTOTIC_SCOPE: finite checks support the lemmas; the infinite limit requires a proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
