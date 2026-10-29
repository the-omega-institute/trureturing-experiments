#!/usr/bin/env python3
"""Exact bounded checks for arXiv:2610.28544v1; no numerical root evaluation."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import comb, lcm, prod


def esymms(x):
    """Independent coefficient recurrence for product_i (1 + x_i t)."""
    c = [F(1)] + [F(0)] * len(x)
    for v in x:
        for j in range(len(x), 0, -1):
            c[j] += v * c[j - 1]
    return c


def sums(x, r):
    return [sum((x[i] for i in I), F(0))
            for I in combinations(range(len(x)), r)]


def root_le(left, right):
    """For nonnegative rational bases, compare A^(1/d) <= B^(1/e)."""
    A, d = left
    B, e = right
    assert A >= 0 and B >= 0 and d > 0 and e > 0
    L = lcm(d, e)
    return A ** (L // d) <= B ** (L // e)


def included(n, left, right):
    """The proposed exhaustive index ranges, including equality of labels."""
    t, a = left
    u, b = right
    if t == 'G' and u == 'G':
        return a >= b
    if t == 'P' and u == 'P':
        return a <= b
    if t == 'G' and u == 'P':
        return a + b >= n + 1
    return a * b <= n


def check_obstructions():
    counts = Counter()
    for n in range(2, 13):
        for a in range(1, n + 1):
            for b in range(1, n + 1):
                if a < b:
                    u = (F(n-a, a) + F(n-b, b)) / 2
                    x = [-u] + [F(1)] * (n-1)
                    c = esymms(x)
                    assert all(c[j] > 0 for j in range(1, a+1))
                    assert c[b] < 0
                    counts['G_not_G'] += 1
                if a > b:
                    u = F(a+b-2, 2)
                    x = [-u] + [F(1)] * (n-1)
                    assert min(sums(x, a)) > 0
                    assert min(sums(x, b)) < 0
                    counts['P_not_P'] += 1
                if a+b <= n:
                    m = n-b
                    eps = F(1, 2 ** (n+1))
                    x = [-eps] * b + [F(1)] * m
                    c = esymms(x)
                    for j in range(1, a+1):
                        terms = [(-eps)**h * comb(b, h) * comb(m, j-h)
                                 for h in range(min(b, j)+1)]
                        assert c[j] == sum(terms)
                        tail = sum(abs(t) for t in terms[1:])
                        assert terms[0] >= 1
                        assert tail <= eps * comb(n, j) <= F(1, 2)
                        assert c[j] >= F(1, 2) > 0
                    assert sum(x[:b]) == -b * eps < 0
                    counts['G_not_P'] += 1
                if a*b > n:
                    u = (F(n-a, a) + (b-1)) / 2
                    x = [-u] + [F(1)] * (n-1)
                    assert min(sums(x, b)) > 0
                    assert esymms(x)[a] < 0
                    counts['P_not_G'] += 1
        for q in [F(-2), F(0), F(1, 7), F(2), F(100, 3)]:
            x = [-q] + [F(1)] * (n-1)
            c = esymms(x)
            for j in range(1, n+1):
                assert c[j] == comb(n-1, j) - q * comb(n-1, j-1)
                assert c[j] == comb(n-1, j-1) * (F(n-j, j) - q)
                counts['one_negative_formula'] += 1
    print('OBSTRUCTION_WITNESSES n=2..12:', dict(sorted(counts.items())))


def check_small_claim_and_source():
    grid = tuple(map(F, [-2, -1, 0, 1, 2])) + (F(-1, 2), F(1, 2))
    grid = tuple(sorted(grid))
    count = Counter()
    for n in range(2, 7):
        labels = [(t, k) for t in ['G', 'P'] for k in range(1, n+1)]
        pairs = [(a, b) for a in labels for b in labels if included(n, a, b)]
        for x in combinations_with_replacement(grid, n):
            count['vectors'] += 1
            c = esymms(x)
            rsums = {k: sums(x, k) for k in range(1, n+1)}
            closed = {}
            mean = {}
            for k in range(1, n+1):
                closed['G', k] = all(c[j] >= 0 for j in range(1, k+1))
                closed['P', k] = min(rsums[k]) >= 0
                mean['G', k] = (c[k] / comb(n, k), k)
                degree = comb(n, k)
                mean['P', k] = (prod(rsums[k]) / F(k)**degree, degree)
            for a, b in pairs:
                if closed[a]:
                    assert closed[b], (n, x, a, b, 'inclusion')
                    assert root_le(mean[a], mean[b]), (n, x, a, b, 'mean')
                    count[a[0]+'_to_'+b[0]] += 1
            # Source Theorem 1.1, (1.6): M_k^2 >= M_{k-1} M_{k+1}.
            for k in range(2, n):
                if closed['P', k-1]:
                    A, d = mean['P', k]
                    B, e = mean['P', k-1]
                    C, f = mean['P', k+1]
                    L = lcm(d, e, f)
                    assert A**(2*L//d) >= B**(L//e) * C**(L//f)
                    count['source_newton'] += 1
            # Source identities Pi_1=sigma_n and Pi_n=sigma_1.
            assert prod(rsums[1]) == c[n]
            assert prod(rsums[n]) == c[1]
            count['source_endpoint_identities'] += 2
    print('SMALL_CASES sorted rational grid {-2,-1,-1/2,0,1/2,1,2}, n=2..6:',
          dict(sorted(count.items())))
    # Concrete fidelity reading: x=(1,2,3), source Theorem 1.1 at k=2.
    # After cubing: (Pi_2/2^3)^2 >= Pi_1 * (Pi_3/3)^3.
    x = (F(1), F(2), F(3))
    lhs = (prod(sums(x, 2))/8)**2
    rhs = prod(sums(x, 1)) * (prod(sums(x, 3))/3)**3
    assert (lhs, rhs) == (F(225, 4), F(48)) and lhs >= rhs
    print(f'FIDELITY Theorem 1.1 n=3 k=2 x=(1,2,3), cubed: {lhs} >= {rhs}')


if __name__ == '__main__':
    check_obstructions()
    check_small_claim_and_source()
    print('ALL_CHECKS_OK')
