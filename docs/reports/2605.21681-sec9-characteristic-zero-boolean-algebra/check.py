#!/usr/bin/env python3
"""Exact finite tests for the proposed Boolean-algebra refutation of arXiv:2605.21681 §9.

Python standard library only. Matrices have integer entries. Primitive integer
row elimination computes rank over Q without approximation. Finite tests do not
prove atomlessness, back-and-forth, invariance under the infinite group, or the
unbounded conclusion; those are explicit obligations of the proposed argument.
"""
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd, comb


def primitive(row):
    g = reduce(gcd, row, 0)
    return [x // g for x in row] if g else row


def rank_q(rows):
    pivots = {}
    for data in rows:
        row = list(data)
        for j, x in enumerate(row):
            if not x:
                continue
            if j in pivots:
                b = pivots[j]
                row = primitive([b[j] * a - x * c for a, c in zip(row, b)])
            else:
                pivots[j] = primitive(row)
                break
    return len(pivots)


def row_masks(n, d):
    return [sum(1 << i for i in s) for j in range(1, min(n, d) + 1)
            for s in combinations(range(n), j)]


def test_matrix(n, d):
    # A tuple gives the same row as its set of distinct coordinates. Repeated
    # coordinates are allowed, so all nonempty subsets of size <= d occur.
    return [[int((t & b) == t) for b in range(1, (1 << n) - 1)]
            for t in row_masks(n, d)]


def fidelity():
    # Source Section 2, Example 5 (vector atoms): in A^2 there is the independent
    # orbit and 2+|F_2| other orbits. GL(3,2) already realizes all five types.
    bases = [(a, b, c) for a, b, c in product(range(1, 8), repeat=3)
             if b != a and c not in (a, b, a ^ b)]
    assert len(bases) == 168
    def action(basis, x):
        r = 0
        for i in range(3):
            if (x >> i) & 1:
                r ^= basis[i]
        return r
    todo = set(product(range(8), repeat=2))
    sizes = []
    while todo:
        x, y = min(todo)
        orb = {(action(b, x), action(b, y)) for b in bases}
        assert orb <= todo
        sizes.append(len(orb))
        todo -= orb
    assert sorted(sizes) == [1, 7, 7, 7, 42]
    print('FIDELITY vector-atoms q=2, dimension=3, GL-order=168: pair-orbits=5, sizes=[1,7,7,7,42]')


def dimensions():
    for k in range(1, 4):
        n = 1 << k
        dim = (1 << n) - 2
        values = []
        for d in range(1, n):
            matrix = test_matrix(n, d)
            r = rank_q(matrix)
            expected = sum(comb(n, j) for j in range(1, d + 1))
            assert r == expected
            values.append(dim - r)
            # Adjacent atom transpositions generate all finite algebra automorphisms.
            rows = set(row_masks(n, d))
            for a in range(n - 1):
                def swap(mask):
                    x, y = (mask >> a) & 1, (mask >> (a + 1)) & 1
                    return mask ^ (((1 << a) | (1 << (a + 1))) if x != y else 0)
                assert {swap(t) for t in rows} == rows
            print(f'DIM_Q k={k} atoms={n} algebra-size={1 << n} V-dim={dim} d={d} rank={r} kernel-dim={dim-r}')
        assert all(x > y for x, y in zip(values, values[1:]))
    print('INVARIANCE: compressed test rows closed under every adjacent atom swap for k=1,2,3')


def refinements():
    checked = 0
    # Exhaust every onto map between atom sets with at most four fine atoms.
    for m in range(1, 5):
        for n in range(1, m + 1):
            for pi in product(range(n), repeat=m):
                if set(pi) != set(range(n)):
                    continue
                for b in range(1, (1 << n) - 1):
                    image = sum(1 << j for j in range(m) if (b >> pi[j]) & 1)
                    assert 0 < image < (1 << m) - 1
                    for d in range(1, 4):
                        for t in product(range(m), repeat=d):
                            lhs = int(all((image >> j) & 1 for j in t))
                            rhs = int(all((b >> pi[j]) & 1 for j in t))
                            assert lhs == rhs
                            checked += 1
                # Pullback injectivity: every coarse tuple has a fine lift.
                for d in range(1, 4):
                    assert {tuple(pi[j] for j in t) for t in product(range(m), repeat=d)} == set(product(range(n), repeat=d))
    print(f'REFINEMENT: all onto atom maps with fine-atoms<=4, d=1,2,3: {checked} basis/tuple equalities; pullbacks injective')


def witnesses():
    for d in range(1, 6):
        n = d + 3
        c, r = 0, n - 1
        coeff = {}
        for s in range(1 << (d + 1)):
            b = 1 << c
            for i in range(d + 1):
                if (s >> i) & 1:
                    b |= 1 << (i + 1)
            coeff[b] = (-1) ** (d + 1 - s.bit_count())
        assert len(coeff) == 1 << (d + 1)
        assert all(0 < b < (1 << n) - 1 for b in coeff)
        def test(t):
            return sum(v for b, v in coeff.items() if all((b >> j) & 1 for j in t))
        vals = [test(t) for t in product(range(n), repeat=d)]
        assert all(v == 0 for v in vals)
        separated = test(tuple(range(1, d + 2)))
        assert separated == 1
        # Integer identities imply these same readings in every characteristic.
        for p in (2, 3, 5, 7):
            assert all(v % p == 0 for v in vals) and separated % p == 1
        print(f'WITNESS d={d} partition-atoms={n} support={len(coeff)} all-{d}-tuple-tests=0 separating-{d+1}-tuple-test=1')
    print('CHARACTERISTICS: integer identities and reductions mod 2,3,5,7 agree')


if __name__ == '__main__':
    fidelity()
    dimensions()
    refinements()
    witnesses()
    print('ALL_CHECKS_OK')
