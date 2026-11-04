#!/usr/bin/env python3
"""Exact finite checks for Li--Liu, arXiv:2609.31982v1, Section 5.4."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial


def connected(n, edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for u, v in edges:
        a, b = find(u), find(v)
        parent[a] = b
    return len({find(i) for i in range(n)}) == 1


def count_words(n, edges):
    """DP records the actual permutation and exactly which edges have appeared.

    Right multiplication by (u v) swaps entries u and v of the permutation.
    A mask is accepted precisely when its spanning support is connected.
    Thus words may omit edges of the ambient graph.
    """
    perms = list(permutations(range(n)))
    index = {p: i for i, p in enumerate(perms)}
    width = 1 << len(edges)
    transitions = []
    for p in perms:
        row = []
        for u, v in edges:
            a = list(p)
            a[u], a[v] = a[v], a[u]
            row.append(index[tuple(a)] * width)
        transitions.append(row)
    dp = {index[tuple(range(n))] * width: 1}
    bits = [1 << i for i in range(len(edges))]
    for _ in range(2 * n - 2):
        out = defaultdict(int)
        for state, count in dp.items():
            p, mask = divmod(state, width)
            for nxt, bit in zip(transitions[p], bits):
                out[nxt + (mask | bit)] += count
        dp = out
    identity_base = index[tuple(range(n))] * width
    total = 0
    for mask in range(width):
        used = [e for i, e in enumerate(edges) if mask & (1 << i)]
        if connected(n, used):
            total += dp.get(identity_base + mask, 0)
    return total, dp.get(identity_base + width - 1, 0)


def weighted_trees(n, edges):
    """Enumerate spanning trees and positive weights, without using classification.

    Orient u<v. Summing nonnegative discrepancies over any initial label interval
    bounds every crossing edge weight by that interval's size, hence by n-1.
    This is a complete finite weight bound, not a numerical cutoff assumption.
    """
    number, total = 0, Fraction(0)
    for tree in combinations(edges, n - 1):
        if not connected(n, tree):
            continue
        for weights in product(range(1, n), repeat=n - 1):
            discrepancy = [1] * n
            numerator = factorial(n)
            for (u, v), w in zip(tree, weights):
                discrepancy[u] -= w
                discrepancy[v] += w
                numerator *= w
            if min(discrepancy) < 0:
                continue
            denominator = 1
            for s in discrepancy:
                denominator *= factorial(s)
            weight = Fraction(numerator, denominator)
            assert weight.denominator == 1
            number += 1
            total += weight
    return number, total


def k11n_edges(n):
    return [(u, v) for u in range(n) for v in range(u + 1, n)
            if v >= n - 2]


def closed_count(n):
    return Fraction(factorial(2*n-2), factorial(n)) * Fraction(
        factorial(2*n-3), factorial(n-2)**2)


def classified_sum(m):
    return sum(comb(m, k) * sum(w * comb(m+2, k+1-w)
               for w in range(1, k+2)) for k in range(m+1))


def main():
    # Fidelity: Section 4.2 explicitly gives 16 factorization trees, weight 96,
    # and N_id(K_4)=2880. These are not inferred from the target conjecture.
    complete4 = list(combinations(range(4), 2))
    tree_number, tree_weight = weighted_trees(4, complete4)
    count, _ = count_words(4, complete4)
    assert (tree_number, tree_weight, count) == (16, 96, 2880)
    print(f'FIDELITY K4 trees={tree_number} weight={tree_weight} N_id={count}')
    # Fidelity: Section 2.2 gives the star count (2n-2)!/n!.
    for n in range(2, 7):
        count, _ = count_words(n, [(i, n-1) for i in range(n-1)])
        expected = Fraction(factorial(2*n-2), factorial(n))
        assert count == expected
        print(f'FIDELITY star n={n} N_id={count}')
    for n in range(3, 7):
        edges = k11n_edges(n)
        count, full_support = count_words(n, edges)
        tree_number, tree_weight = weighted_trees(n, edges)
        weight_formula = Fraction(factorial(2*n-3), factorial(n-2)**2)
        assert count == closed_count(n)
        assert tree_weight == classified_sum(n-2) == weight_formula
        assert count == Fraction(factorial(2*n-2), factorial(n)) * tree_weight
        print(f'TARGET n={n} N_id={count} full_support={full_support} '
              f'trees={tree_number} weight={tree_weight}')
        if n == 5:
            assert full_support == 1920  # Shapiro, arXiv:2610.14373v1, Theorem 6.
            assert count == 47040 and full_support < count
    # Scalar identities alone do not settle the transposition counting claim.
    for m in range(101):
        N = m + 1
        W = classified_sum(m)
        assert W == sum(w * comb(2*N, N-w) for w in range(1, N+1))
        assert W == N * comb(2*N-1, N-1)
        assert W == Fraction(factorial(2*m+1), factorial(m)**2)
        for w in range(1, N+1):
            r = N-w
            previous = comb(2*N-1, r-1) if r > 0 else 0
            assert w * comb(2*N, r) == N * (comb(2*N-1, r)-previous)
    print('SCALAR identities: m=0..100; every telescoping summand checked exactly')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
