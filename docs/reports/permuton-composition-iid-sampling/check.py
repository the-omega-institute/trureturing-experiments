#!/usr/bin/env python3
"""Exact finite checks for Gorski's Section 11 iid-sampling question.

All evidential quantities are integers or fractions.Fraction. No asymptotic
claim is inferred from the finite checks. Run with Python 3; no dependencies.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial


def compose(sigma, tau):
    return tuple(sigma[tau[i]] for i in range(len(tau)))


def matrix(sigma):
    n = len(sigma)
    return tuple(tuple(F(int(sigma[i] == j), n) for j in range(n)) for i in range(n))


def markov_product(p, q):
    """Cell masses: Q acts first, hence n * Q_matrix * P_matrix."""
    n = len(p)
    return tuple(tuple(n * sum(q[i][k] * p[k][j] for k in range(n))
                       for j in range(n)) for i in range(n))


def induced(points):
    n = len(points)
    x_order = sorted(range(n), key=lambda i: (points[i][0], i))
    y_order = sorted(range(n), key=lambda i: (points[i][1], i))
    y_rank = {label: rank for rank, label in enumerate(y_order)}
    return tuple(y_rank[label] for label in x_order)


def match_middle(p, q):
    p_order = sorted(p, key=lambda z: z[0])
    q_order = sorted(q, key=lambda z: z[1])
    return tuple((a[0], b[1]) for a, b in zip(q_order, p_order))


def fidelity():
    """Section 3: checkerboards multiply as permutations of the same size.

    Section 1 also states that the product of independent uniform permutations
    is uniform; counting all pairs tests this finite distributional identity.
    """
    count = 0
    for n in range(1, 6):
        ps = tuple(permutations(range(n)))
        matrices = {s: matrix(s) for s in ps}
        products = Counter()
        u = tuple(tuple(F(1, n * n) for _ in range(n)) for _ in range(n))
        for s, t in product(ps, repeat=2):
            st = compose(s, t)
            assert markov_product(matrices[s], matrices[t]) == matrices[st]
            products[st] += 1
            count += 1
        assert len(products) == factorial(n)
        assert set(products.values()) == {factorial(n)}
        for s in ps:
            assert markov_product(matrices[s], u) == u
            assert markov_product(u, matrices[s]) == u
    print(f"FIDELITY: checkerboard composition, uniform-product law, absorption; {count} pairs, 1 <= n <= 5")


def ranking_identity():
    """Every coordinate-order type and every sample relabeling for n <= 3."""
    count = 0
    for n in range(1, 4):
        ps = tuple(permutations(range(n)))
        for s, t, lp, lq in product(ps, repeat=4):
            # Different middle arrays, with matching ranks but different values.
            p0 = tuple((F(5 * i + 1, 5 * n), F(5 * s[i] + 2, 5 * n))
                       for i in range(n))
            q0 = tuple((F(5 * i + 3, 5 * n), F(5 * t[i] + 4, 5 * n))
                       for i in range(n))
            p = tuple(p0[i] for i in lp)
            q = tuple(q0[i] for i in lq)
            assert induced(p) == s and induced(q) == t
            assert induced(match_middle(p, q)) == compose(s, t)
            count += 1
    print(f"SAMPLING IDENTITY: all {count} labeled order types, 1 <= n <= 3")


def order_cdf(n, k, t):
    """P(U_(k) <= t) = P(Binomial(n,t) >= k), k is one-based."""
    return sum(F(comb(n, j)) * t**j * (1 - t)**(n - j)
               for j in range(k, n + 1))


def beta_cdf(n, k, t):
    """Integrate the exact order-statistic density by binomial expansion."""
    coefficient = F(factorial(n), factorial(k - 1) * factorial(n - k))
    return coefficient * sum(F((-1)**j * comb(n - k, j), k + j) * t**(k + j)
                             for j in range(n - k + 1))


def matched_step_expectation(n, f, g):
    """Independent middle arrays: E[n^-1 sum f(U_(k))g(V_(k))]."""
    m = len(f)
    assert len(g) == m
    weights = []
    for k in range(1, n + 1):
        cdf = [order_cdf(n, k, F(j, m)) for j in range(m + 1)]
        for j, value in enumerate(cdf):
            assert value == beta_cdf(n, k, F(j, m))
        row = [cdf[j + 1] - cdf[j] for j in range(m)]
        assert sum(row) == 1 and all(w >= 0 for w in row)
        weights.append(row)
    for j in range(m):
        assert sum(row[j] for row in weights) == F(n, m)
    return sum(sum(a * w for a, w in zip(f, row)) *
               sum(b * w for b, w in zip(g, row))
               for row in weights) / n


def sampling_expectations():
    sizes = (1, 2, 4, 8, 16, 32, 64)
    # Singular diagonal and antidiagonal permutons: lower-half conditional
    # probabilities are respectively (1,0) and (0,1).
    for n in sizes:
        same = matched_step_expectation(n, (F(1), F(0)), (F(1), F(0)))
        opposite = matched_step_expectation(n, (F(1), F(0)), (F(0), F(1)))
        # Independent Binomial(n,1/2) threshold counts give the same answer.
        pmf = [F(comb(n, j), 2**n) for j in range(n + 1)]
        independent_count_formula = sum(F(min(a, b), n) * pmf[a] * pmf[b]
                                        for a, b in product(range(n + 1), repeat=2))
        assert same == independent_count_formula
        assert same + opposite == F(1, 2)
        error = F(1, 2) - same
        assert error >= 0 and error**2 <= F(1, 8 * n)
    print("SINGULAR FAMILY: Delta/Delta and Delta/antidiagonal exact expectations; n = 1,2,4,8,16,32,64; error^2 <= 1/(8n)")
    # Noncommuting rational stationary kernels on three equal cells. Each row
    # sums to 1, each column sums to 1. P's output event is column 0, Q's input
    # event is row 0. The continuum rectangle mass is (Q*P)[0,0]/3.
    p = ((F(1, 2), F(1, 2), F(0)),
         (F(0), F(1, 2), F(1, 2)),
         (F(1, 2), F(0), F(1, 2)))
    q = ((F(1, 3), F(2, 3), F(0)),
         (F(2, 3), F(0), F(1, 3)),
         (F(0), F(1, 3), F(2, 3)))
    for a in (p, q):
        assert all(sum(row) == 1 for row in a)
        assert all(sum(a[i][j] for i in range(3)) == 1 for j in range(3))
    f, g = tuple(p[i][0] for i in range(3)), q[0]
    target = sum(a * b for a, b in zip(f, g)) / 3
    assert target == F(1, 18)
    assert sum(p[0][i] * q[i][0] for i in range(3)) / 3 == F(1, 6)
    for n in sizes:
        expectation = matched_step_expectation(n, f, g)
        # A finite check of an O(n^-1/2) envelope, not a proof of the limit.
        assert n * (expectation - target)**2 <= 1
    print("NONCOMMUTING STEP KERNELS: Q-then-P rectangle = 1/18, reversed = 1/6; exact order-statistic expectations checked")


def main():
    fidelity()
    ranking_identity()
    sampling_expectations()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
