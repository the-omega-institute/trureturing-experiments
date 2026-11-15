#!/usr/bin/env python3
"""Exact finite certificates for the uniform finite GW/fGW local model.

All arithmetic used as evidence is integer or Fraction arithmetic. The solver
certifies the global optimum over the continuous Birkhoff polytope, rather than
asserting that testing its vertices alone solves a quadratic program. A checked
lower bound on every ordered permutation cross-cost, together with Birkhoff's
theorem, certifies all real convex combinations simultaneously. This is a finite
check of the proposed local-equality lemma, not a consistency proof.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial, isqrt


def metric(n, values):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for (i, j), v in zip(combinations(range(n), 2), values):
        a[i][j] = a[j][i] = F(v)
    assert all(a[i][j] <= a[i][k] + a[k][j]
               for i, j, k in product(range(n), repeat=3))
    return a


def rational_power(x, p, q):
    if p % q == 0:
        return x ** (p // q)
    assert (p, q) == (1, 2)
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    assert a*a == x.numerator and b*b == x.denominator
    return F(a, b)


def kernel(a, b, u, v, p, q, alpha):
    n = len(a)
    return [[[[((abs(a[i][k] - b[j][l]) ** p) if alpha is None else
                rational_power(alpha * abs(a[i][k] - b[j][l]) ** q
                 + (1-alpha) * abs(u[i] - v[j]) ** q, p, q))
               for l in range(n)] for k in range(n)]
             for j in range(n)] for i in range(n)]


def cross_table(c, perms):
    n = len(c)
    return [[sum((c[i][s[i]][k][t[k]] for i in range(n) for k in range(n)), F(0))
             / (n*n) for t in perms] for s in perms]


def objective(c, weights, perms):
    n = len(c)
    P = [[sum((w for w, s in zip(weights, perms) if s[i] == j), F(0))
          for j in range(n)] for i in range(n)]
    assert all(sum(row) == 1 for row in P)
    assert all(sum(P[i][j] for i in range(n)) == 1 for j in range(n))
    return sum((c[i][j][k][l]*P[i][j]*P[k][l]
                for i, j, k, l in product(range(n), repeat=4)), F(0))/(n*n)


def main():
    cases = cross_checks = mixture_checks = fidelity_checks = 0
    modes = [(p, 1, None) for p in range(1, 5)] + [
        (1, 1, F(1, 2)), (2, 1, F(1, 2)),
        (2, 2, F(1, 2)), (4, 2, F(1, 2)), (1, 2, F(1, 4))]
    for n in range(1, 5):
        perms = list(permutations(range(n)))
        per_n = 0
        for edges in product((1, 2), repeat=n*(n-1)//2):
            a = metric(n, edges)
            eps = F(1, 100*n*n)
            b = [[a[i][j] + (eps if i != j else 0) for j in range(n)]
                 for i in range(n)]
            u = [F(i, n) for i in range(n)]
            v = [x + eps for x in u]
            for p, q, alpha in modes:
                features = ([F(0)]*n, [F(0)]*n) if p < q else (u, v)
                c = kernel(a, b, *features, p, q, alpha)
                table = cross_table(c, perms)
                m = min(table[i][i] for i in range(len(perms)))
                if n == 1:
                    optimum = m
                else:
                    delta = min(b[i][j] for i in range(n) for j in range(n) if i != j)
                    gap = delta**p / (n*n)
                    if alpha is not None:
                        gap *= rational_power(alpha, p, q)
                    assert m < gap
                    for i, s in enumerate(perms):
                        for j, t in enumerate(perms):
                            if i != j:
                                assert table[i][j] >= gap
                                # The diagonal source-pair term alone gives the gap.
                                witness = next(k for k in range(n) if s[k] != t[k])
                                assert c[witness][s[witness]][witness][t[witness]] / (n*n) >= gap
                                cross_checks += 1
                    # Every entry is >= m. For every probability vector lambda,
                    # lambda^T table lambda >= m*(sum lambda)^2 = m.
                    assert all(z >= m for row in table for z in row)
                    optimum = m  # attained by a minimizing permutation
                # Reproduce the source's GW_2 <= permutation quotient inequality.
                if alpha is None and p == 2:
                    assert optimum <= min(table[i][i] for i in range(len(perms)))
                    # Source's nonnegative cross-term / n! lower bound, squared.
                    assert all(z >= 0 for row in table for z in row)
                    assert optimum >= min(table[i][i] for i in range(len(perms)))/factorial(n)
                    fidelity_checks += 1
                # Independent summation of the continuous objective for rational
                # points, compared to the Birkhoff cross-term expansion.
                vectors = [[F(1, len(perms))] * len(perms)]
                if len(perms) > 1:
                    for weight in (F(1, 4), F(1, 2), F(3, 4)):
                        z = [F(0)] * len(perms)
                        z[0], z[-1] = weight, 1-weight
                        vectors.append(z)
                for weights in vectors:
                    lhs = objective(c, weights, perms)
                    rhs = sum((weights[i]*weights[j]*table[i][j]
                               for i, j in product(range(len(perms)), repeat=2)), F(0))
                    assert lhs == rhs and lhs >= optimum
                    mixture_checks += 1
                cases += 1
                per_n += 1
        print(f'n={n}: exact globally certified local GW/fGW instances={per_n}')
    # A separate exact complete solution for every two-point source metric in
    # this integer box. Q(t)=m+2*a*b*t*(1-t), 0<=t<=1, for GW_2.
    for a, b in product(range(1, 6), repeat=2):
        A, B = metric(2, (a,)), metric(2, (b,))
        T = cross_table(kernel(A, B, [0, 0], [0, 0], 2, 1, None), [(0, 1), (1, 0)])
        m = F((a-b)**2, 2)
        assert T[0][0] == T[1][1] == m
        assert T[0][1] == T[1][0] == F(a*a+b*b, 2)
        assert T[0][1]-m == a*b > 0
        assert (m == 0) == (a == b)  # source zero iff measure-preserving isometry
        fidelity_checks += 1
    print(f'local_certificates={cases}; ordered_cross_certificates={cross_checks}')
    print(f'independent_objective_expansions={mixture_checks}; source_fidelity_checks={fidelity_checks}')
    print('Scope: GW p=1,2,3,4; fused (1,1),(2,1),(2,2),(4,2) at alpha=1/2; (1,2) at alpha=1/4 with zero features; d=1.')
    print('Continuum certificate uses Birkhoff; no classifier limit or arbitrary-real exponent is numerically certified.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
