#!/usr/bin/env python3
"""Exact independent checks of Combe's model and Chapoton's A277506 conjecture.

The cover relation is derived from all componentwise intervals, without using
our proposed cover classification or continuation recurrence. Python 3 only;
all arithmetic is integer or fractions.Fraction. This is finite evidence.
"""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import factorial

OEIS = [1, 2, 8, 51, 460, 5425, 79206, 1377985, 27801096,
        637630353, 16376303530, 465451009441, 14501512561548]


def triwords(d):
    return [w for w in product(range(3), repeat=d)
            if w[0] != 2 and all(not (w[i] == 0 and w[j] == 1)
                                  for i in range(d) for j in range(i + 1, d))]


def componentwise_covers(words):
    size = len(words)
    above, below = [0] * size, [0] * size
    for i, u in enumerate(words):
        for j, v in enumerate(words):
            if i != j and all(a <= b for a, b in zip(u, v)):
                above[i] |= 1 << j
                below[j] |= 1 << i
    # Intersection is exactly the set of strict intermediate words.
    return [[j for j in range(size)
             if above[i] & (1 << j) and above[i] & below[j] == 0]
            for i in range(size)]


def proposed_covers(w):
    result = set()
    first_zero = next((i for i, a in enumerate(w) if a == 0), None)
    for i, a in enumerate(w):
        v = list(w)
        if a == 0:
            v[i] = 1 if i == first_zero else 2
        elif a == 1 and i != 0:
            v[i] = 2
        else:
            continue
        result.add(tuple(v))
    return result


@cache
def F(z, k):
    if z == 0:
        return factorial(k)
    return F(z - 1, k + 1) + (z - 1) * F(z - 1, k) + (k * F(z, k - 1) if k else 0)


def exp_series(f, degree):
    assert f[0] == 0
    e = [Q(1)] + [Q(0)] * degree
    for n in range(1, degree + 1):
        e[n] = sum((k * f[k] * e[n - k] for k in range(1, n + 1)), Q(0)) / n
    return e


def multiply(f, g, degree):
    return [sum((f[k] * g[n - k] for k in range(n + 1)), Q(0))
            for n in range(degree + 1)]


def inverse_unit(f, degree):
    assert f[0] == 1
    g = [Q(1)] + [Q(0)] * degree
    for n in range(1, degree + 1):
        g[n] = -sum((f[k] * g[n - k] for k in range(1, n + 1)), Q(0))
    return g


def main():
    examples = {1: {'0', '1'}, 2: {'00', '02', '10', '11', '12'},
                3: {'000', '020', '002', '100', '022', '110', '102', '120',
                    '111', '121', '112', '122'}}
    for d, expected in examples.items():
        assert {''.join(map(str, w)) for w in triwords(d)} == expected
    print('FIDELITY: exact triword sets in Combe Section 1.1 (journal Section 2.1), d=1,2,3: OK')
    totals = []
    for d in range(1, 9):
        words = triwords(d)
        assert 4 * len(words) == 2 ** d * (d + 3)
        covers = componentwise_covers(words)
        top = words.index((1,) + (2,) * (d - 1))
        bottom = words.index((0,) * d)
        counts, lengths = [0] * len(words), [set() for _ in words]
        counts[top], lengths[top] = 1, {0}
        for i in sorted(range(len(words)), key=lambda i: sum(words[i]), reverse=True):
            if i != top:
                counts[i] = sum(counts[j] for j in covers[i])
                lengths[i] = {ell + 1 for j in covers[i] for ell in lengths[j]}
            assert {words[j] for j in covers[i]} == proposed_covers(words[i])
            if words[i][0] == 1:
                assert counts[i] == F(words[i].count(0), words[i][1:].count(1))
        assert max(lengths[bottom]) == 2 * d - 1  # Combe's longest-chain lemma.
        assert counts[bottom] == OEIS[d - 1]
        totals.append(counts[bottom])
        print(f'd={d}: vertices={len(words)}, maximal_chains={counts[bottom]}, '
              f'chain_edge_lengths={sorted(lengths[bottom])}; covers/state recurrence/longest length: OK')
        if d == 2:
            assert lengths[bottom] == {2, 3} and counts[bottom] == 2
    print('FIDELITY: 4*#Tr(d)=2^d*(d+3), and longest chain length 2d-1, d=1..8: OK')
    print('SETTLEMENT: interval-derived maximal-chain counts, n=0..7:', totals)
    degree = len(OEIS) - 1
    T = [Q(0)] * (degree + 1)
    for n in range(1, degree + 1):
        T[n] = exp_series(T, degree)[n - 1]
    e = exp_series(T, degree)
    assert T == [Q(0)] + e[:-1]  # T = x exp(T) mod x^(degree+1).
    # W(x) = -T(-x), and W(x) exp(W(x)) = x.
    W = [-((-1) ** n) * T[n] for n in range(degree + 1)]
    assert multiply(W, exp_series(W, degree), degree) == [Q(0), Q(1)] + [Q(0)] * (degree - 1)
    one_minus_T = [-v for v in T]
    one_minus_T[0] += 1
    one_minus_x = [Q(1), Q(-1)] + [Q(0)] * (degree - 1)
    B = inverse_unit(one_minus_T, degree)
    C = inverse_unit(multiply(one_minus_T, one_minus_x, degree), degree)
    for n, expected in enumerate(OEIS):
        assert factorial(n) * C[n] == expected
        assert factorial(n) * B[n] == F(n, 0)
        if n:
            assert expected == n * OEIS[n - 1] + n ** n
            assert F(n, 0) == n ** n
    print(f'FORMAL SERIES (exact rational truncation): W exp(W)=x and defining e.g.f. through x^{degree}: OK')
    print(f'OEIS recurrence and F(n,0) versus series, n=0..{degree}: OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
