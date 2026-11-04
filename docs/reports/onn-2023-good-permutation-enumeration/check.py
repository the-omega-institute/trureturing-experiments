#!/usr/bin/env python3
"""Exact independent checks for Onn's unnumbered referee question.

Only Python's standard library and integer arithmetic are used. Enumeration
assigns every rank once and prunes only after an already completed square fails
the literal source definition. Rejected prefixes account for all their factorial
many completions. No claimed count or chain condition is used to prune.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
from math import comb, factorial


def good_square(a, b, c, d):
    return max(a, b, c, d) in (a, d) and min(a, b, c, d) in (b, c)


def good(a, m, n):
    return all(good_square(a[i*n+j], a[i*n+j+1],
                           a[(i+1)*n+j], a[(i+1)*n+j+1])
               for i in range(m-1) for j in range(n-1))


def transpose(a, m, n):
    return tuple(a[i*n+j] for j in range(n) for i in range(m))


def inc(chain):
    return all(a < b for a, b in zip(chain, chain[1:]))


def chain_splits(a, n):
    top, bottom = a[:n], a[n:]
    return [k for k in range(n+1)
            if inc(bottom[:k] + top[:k][::-1])
            and inc(top[k:][::-1] + bottom[k:])]


def types33(a):
    aa, b, c, d, e, f, g, h, i = a
    return (g < d < e < b < aa and c < f < e < h < i,
            c < b < e < d < aa and g < h < e < f < i)


def proposed(m, n):
    s, t = sorted((m, n))
    if s == 1:
        return factorial(t)
    if s == 2:
        return 2 ** (2*t-1)
    return 72 if s == t == 3 else 0


def empty_square_count(size):
    # Exhaustive permutation search with identical suffix subproblems memoized.
    # The recurrence sums every possible next unused label, without a count formula.
    full = (1 << size) - 1

    @lru_cache(None)
    def visit(used):
        if used == full:
            return 1
        return sum(visit(used | (1 << k)) for k in range(size)
                   if not used & (1 << k))
    return visit(0), visit.cache_info().currsize


def enumerate_good(m, n):
    size = m*n
    if min(m, n) == 1:
        count, states = empty_square_count(size)
        return count, None, Counter(), states, 0
    # Complete adjacent squares as early as possible in either orientation.
    slots = ([i*n+j for j in range(n) for i in range(m)] if m <= n
             else list(range(size)))
    squares = [(i*n+j, i*n+j+1, (i+1)*n+j, (i+1)*n+j+1)
               for i in range(m-1) for j in range(n-1)]
    position = {p: k for k, p in enumerate(slots)}
    closing = [[] for _ in slots]
    for square in squares:
        closing[max(position[p] for p in square)].append(square)
    values = [0]*size
    accepted = set()
    splits = Counter()
    nodes = rejected_completions = 0

    def visit(depth, available):
        nonlocal nodes, rejected_completions
        nodes += 1
        if depth == size:
            a = tuple(values)
            assert good(a, m, n)
            assert good(transpose(a, m, n), n, m)
            if m == 2 or n == 2:
                oriented = a if m == 2 else transpose(a, m, n)
                ks = chain_splits(oriented, max(m, n))
                assert len(ks) == 1, (m, n, a, ks)
                splits[ks[0]] += 1
            if m == n == 3:
                ty = types33(a)
                assert sum(ty) == 1 and a[4] == 5
                splits['A' if ty[0] else 'B'] += 1
            accepted.add(a)
            return
        p = slots[depth]
        remaining = available
        while remaining:
            bit = remaining & -remaining
            remaining -= bit
            values[p] = bit.bit_length()
            if all(good_square(*(values[q] for q in sq)) for sq in closing[depth]):
                visit(depth+1, available ^ bit)
            else:
                rejected_completions += factorial(size-depth-1)
        values[p] = 0

    visit(0, (1 << size)-1)
    assert len(accepted) + rejected_completions == factorial(size)
    return len(accepted), accepted, splits, nodes, rejected_completions


def generate_chains(n):
    arrays = set()
    for k in range(n+1):
        count = 0
        for chosen in combinations(range(1, 2*n+1), 2*k):
            chosen = tuple(chosen)
            rest = tuple(x for x in range(1, 2*n+1) if x not in chosen)
            top = chosen[k:][::-1] + rest[:n-k][::-1]
            bottom = chosen[:k] + rest[n-k:]
            a = top + bottom
            assert chain_splits(a, n) == [k]
            assert good(a, 2, n)
            assert a not in arrays
            arrays.add(a)
            count += 1
        assert count == comb(2*n, 2*k)
    return arrays


def main():
    # Fidelity to Theorem 1.2(2), (3) and its displayed additional examples.
    for n in range(1, 7):
        sigma = tuple(range(n, 0, -1)) + tuple(range(n+1, 2*n+1))
        assert good(sigma, 2, n)
    for sigma, m, n in [((8,7,1,4,5,3,2,6,9),3,3),
                         ((9,6,2,3,5,4,1,7,8),3,3),
                         ((4,3,1,2),2,2)]:
        assert good(sigma, m, n)
    # A rank-universality fidelity check: each distinct nondecreasing nine-tuple
    # in {0,1,2} yields a supmodular A^sigma in Theorem 1.2(3).
    sigma = (8,7,1,4,5,3,2,6,9)
    monotone_inputs = 0
    for seq in combinations_with_replacement(range(3), 9):
        a = [seq[k-1] for k in sigma]
        assert all(a[i*3+j]+a[r*3+s] >= a[i*3+s]+a[r*3+j]
                   for i in range(3) for r in range(i+1,3)
                   for j in range(3) for s in range(j+1,3))
        monotone_inputs += 1
    print(f'FIDELITY Theorem 1.2 examples: OK; ranked supmodular inputs={monotone_inputs}', flush=True)
    # Independent complete 2x2 local truth table.
    for a,b,c,d in permutations(range(1,5)):
        signs = (a>c, b>d)
        predicted = ((c<d<b<a) if signs==(True,True) else
                     (b<a<c<d) if signs==(False,False) else
                     signs==(True,False))
        assert good_square(a,b,c,d) == predicted
    print('LOCAL_CASES all 24 rank arrangements: OK', flush=True)
    assert generate_chains(1) == {(1, 2), (2, 1)}
    results = {}
    for m in range(1,13):
        for n in range(1,12//m+1):
            count, arrays, by_split, nodes, rejected = enumerate_good(m,n)
            assert count == proposed(m,n), (m,n,count,proposed(m,n))
            if min(m,n)==2:
                t=max(m,n)
                generated=generate_chains(t)
                normalized=arrays if m==2 else {transpose(a,m,n) for a in arrays}
                assert normalized == generated
                assert by_split == Counter({k:comb(2*t,2*k) for k in range(t+1)})
            if m==n==3:
                assert by_split == Counter(A=36,B=36)
                # Independently classify all 9! complete arrays without pruning.
                total=0
                for a in permutations(range(1,10)):
                    ty=types33(a)
                    assert good(a,3,3) == any(ty)
                    total += any(ty)
                assert total==72
            results[m,n]=count
            print(f'ENUM {m}x{n}: count={count}; visited={nodes}; rejected_completions={rejected}; classes={dict(by_split)}', flush=True)
    assert set(results) == {(m,n) for m in range(1,13) for n in range(1,13) if m*n<=12}
    assert all(results[m,n]==results[n,m] for m,n in results)
    print(f'ALL_SMALL_CASES ordered_dimensions={len(results)}, m*n<=12: OK', flush=True)
    print('TWO_CHAINS exact set equality in both orientations for widths 1..6: OK', flush=True)
    print('3x3 complete 9! classification: A=36, B=36, center=5: OK', flush=True)
    print('3x4 and 4x3 count=0, agreeing with source Theorem 1.3: OK', flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
