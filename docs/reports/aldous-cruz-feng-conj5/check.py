#!/usr/bin/env python3
"""Exact checks for the preregistered attracting rank-cycle family.

Only integers and fractions are used. The uniform family certificate checks
the arithmetic and algebra of the analytical argument; Hoeffding's inequality
and the contraction theorem are mathematical inputs, not numerically tested
substitutes for a proof. This program is not a Lean proof.
"""

from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def binomial_tail(k, j, p):
    """P(Bin(k,p) >= j), computed over one integer denominator."""
    a, b = p.numerator, p.denominator
    return Q(sum(comb(k, r) * a**r * (b-a)**(k-r)
                 for r in range(j, k+1)), b**k)


def binomial_atom(k, r, p):
    a, b = p.numerator, p.denominator
    return Q(comb(k, r) * a**r * (b-a)**(k-r), b**k)


def rank_kernel(orders, theta, j, k):
    n = len(theta)
    rows = []
    for order in orders:
        row = [Q(0)] * n
        cumulative, previous = Q(0), Q(0)
        for y in order:
            cumulative += theta[y]
            current = binomial_tail(k, j, cumulative)
            row[y] = current - previous
            previous = current
        assert sum(row) == 1 and all(x >= 0 for x in row)
        rows.append(row)
    return rows


def vecmul(v, a):
    return [sum(v[i]*a[i][j] for i in range(len(v)))
            for j in range(len(a[0]))]


def matmul(a, b):
    return [vecmul(row, b) for row in a]


def madd(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def stationary(a):
    """Independent exact linear solve, with normalization as the last equation."""
    n = len(a)
    system = [[a[i][j] - Q(i == j) for i in range(n)] + [Q(0)]
              for j in range(n-1)]
    system.append([Q(1)]*n + [Q(1)])
    for c in range(n):
        pivot = next(r for r in range(c, n) if system[r][c])
        system[c], system[pivot] = system[pivot], system[c]
        v = system[c][c]
        system[c] = [x/v for x in system[c]]
        for r in range(n):
            if r != c:
                v = system[r][c]
                system[r] = [x-v*y for x, y in zip(system[r], system[c])]
    p = [row[-1] for row in system]
    assert sum(p) == 1 and min(p) >= 0 and vecmul(p, a) == p
    return p


def fidelity():
    # arXiv:2403.18153v1, section 5.1, equation (R2) and displayed exact K.
    ranks = [[1, 2, 4, 3], [2, 1, 3, 4],
             [4, 2, 1, 3], [2, 4, 3, 1]]
    orders = [sorted(range(4), key=lambda y: row[y]) for row in ranks]
    theta = [Q(1, 6), Q(1, 6), Q(2, 6), Q(2, 6)]
    expected = [[Q(11, 36), Q(1, 4), Q(1, 9), Q(1, 3)],
                [Q(1, 4), Q(11, 36), Q(1, 3), Q(1, 9)],
                [Q(1, 36), Q(7, 36), Q(5, 9), Q(2, 9)],
                [Q(7, 36), Q(1, 36), Q(2, 9), Q(5, 9)]]
    assert rank_kernel(orders, theta, 1, 2) == expected
    assert vecmul(theta, expected) == theta
    assert stationary(expected) == theta
    assert stationary(rank_kernel(orders, theta, 2, 2)) == theta
    print('FIDELITY: source R2 transition matrix and invariant law reproduced exactly')


def sampling_crosscheck():
    cases = 0
    for points in [(0, 1), (0, 1, 3), (0, 4, 7, 12)]:
        n = len(points)
        orders = [sorted(range(n), key=lambda y: abs(points[x]-points[y]))
                  for x in range(n)]
        distributions = [[Q(1, n)]*n,
                         [Q(i+1, n*(n+1)//2) for i in range(n)],
                         [Q(1)] + [Q(0)]*(n-1)]
        for k in range(2, 6):
            for theta in distributions:
                direct = [[[Q(0)]*n for _ in range(n)] for _ in range(k)]
                for samples in product(range(n), repeat=k):
                    mass = Q(1)
                    for y in samples:
                        mass *= theta[y]
                    for x in range(n):
                        # Deterministic sample-index tie breaking has the same
                        # law as uniform tie breaking for iid samples. Repeated
                        # equal-distance samples here are the same state.
                        indices = sorted(range(k), key=lambda i:
                                         (abs(points[x]-points[samples[i]]), i))
                        for r, i in enumerate(indices):
                            direct[r][x][samples[i]] += mass
                for j in range(1, k+1):
                    assert direct[j-1] == rank_kernel(orders, theta, j, k)
                    cases += 1
    print(f'SAMPLING: tuple enumeration equals binomial-tail kernel in {cases} exact cases')


def cone_and_cycle():
    # Positive u=a-b, v=a+b-c, w=c-a parameterize the entire open cone.
    # Then b=v+w, a=u+v+w, c=u+v+2w.
    points = [(0, 0, 0), (1, 1, 1), (1, 2, 2), (2, 3, 4)]
    subtract = lambda x, y: tuple(a-b for a, b in zip(x, y))
    positive = lambda x: all(a >= 0 for a in x) and any(a > 0 for a in x)
    distances = [[subtract(points[max(x, y)], points[min(x, y)])
                  for y in range(4)] for x in range(4)]
    orders = [[0, 1, 2, 3], [1, 2, 0, 3], [2, 1, 3, 0], [3, 2, 1, 0]]
    for row, order in zip(distances, orders):
        for left, right in zip(order, order[1:]):
            assert positive(subtract(row[right], row[left]))
    # b < a < c < a+b < b+c < a+b+c.
    all_distances = [distances[x][y]
                     for x, y in [(1, 2), (0, 1), (2, 3), (0, 2), (1, 3), (0, 3)]]
    for left, right in zip(all_distances, all_distances[1:]):
        assert positive(subtract(right, left))
    concrete = [0, 4, 7, 12]
    assert sorted(abs(concrete[x]-concrete[y])
                  for x in range(4) for y in range(x+1, 4)) == [3, 4, 5, 7, 8, 12]
    assert [row[2] for row in orders] == [2, 0, 3, 1]
    c = [[Q(y == orders[x][2]) for y in range(4)] for x in range(4)]
    identity = [[Q(i == j) for j in range(4)] for i in range(4)]
    pi = [[Q(1, 4)]*4 for _ in range(4)]
    c2, c3 = matmul(c, c), matmul(matmul(c, c), c)
    assert matmul(c3, c) == identity
    assert madd(madd(identity, c), madd(c2, c3)) == scale(4, pi)
    r = scale(Q(1, 4), madd(madd(scale(3, identity), scale(2, c)), c2))
    assert matmul(madd(identity, scale(-1, c)), r) == madd(identity, scale(-1, pi))
    assert all(x >= 0 for row in r for x in row)
    assert [sum(r[i][j] for i in range(4)) for j in range(4)] == [Q(3, 2)]*4
    print('FAMILY: symbolic positive-coefficient certificate covers 0<b<a<c<a+b')
    print('CYCLE: 0->7->12->4->0; (I-C)R=I-Pi and column norm(R)=3/2')
    return orders, c, r


def uniform_arithmetic():
    # The analytical binomial Hoeffding bound is exp(-2*N*gap^2).
    assert Q(2)*(Q(5, 8)-Q(17, 32))**2 == Q(9, 512)
    assert Q(2)*(Q(47, 64)-Q(5, 8))**2 == Q(49, 2048)
    assert Q(49, 2048) >= Q(9, 512)
    # t(k)=(j-1)/(k-1)=5/8-3/(8(k-1)); increasing to 5/8.
    t_min = Q(5, 8) - Q(3, 8*48)
    assert t_min-Q(17, 32) == Q(11, 128)
    assert t_min-Q(17, 64) >= Q(11, 128)
    assert Q(47, 64)-Q(5, 8) >= Q(11, 128)
    assert 2*Q(11, 128)**2 == Q(121, 8192)
    # F' = k*P(Bin(k-1,p)=j-1). The maximum adjacent-prefix
    # directional coefficient on the sum-zero hyperplane is 3.
    assert max(1, 1+2, 2+1, 1) == 3
    k0 = 1024
    assert Q(9*k0, 512) == 18
    assert Q(121*(k0-1), 8192) > 15
    assert 1-Q(121*k0, 8192) < 0
    # The last inequality certifies the sign of d/dx[x*exp(-121(x-1)/8192)]
    # on x>=1024: its remaining factor is positive. e > 2 follows
    # already from the first three terms of its positive power series.
    assert sum(Q(1, factorial(i)) for i in range(3)) > 2
    eps_upper = Q(2, 2**18)
    beta_upper = Q(k0, 2**15)
    lips_kernel_upper = 3*beta_upper
    assert eps_upper == Q(1, 2**17)
    assert lips_kernel_upper == Q(3, 32)
    assert 1-6*eps_upper > Q(3, 4)
    assert Q(3, 2)*eps_upper < Q(1, 64)
    assert (Q(3, 2)*lips_kernel_upper)/Q(3, 4) == Q(3, 16)
    assert (Q(3, 2)*lips_kernel_upper)/(1-6*eps_upper) < Q(3, 16)
    # A coordinate-chart cube of radius 1/256 lies in B: the fourth
    # coordinate deviates by at most 3/256 < 1/64.
    assert 3*Q(1, 256) < Q(1, 64)
    cube_volume = (2*Q(1, 256))**3
    assert cube_volume == Q(1, 2**21) > 0
    print('UNIFORM: k=8m, m>=128: eps<1/131072, L<3/32, stationary Lipschitz<3/16')
    print('BASIN: self-map radius<1/64; fixed-point coordinates>=15/64; chart volume>=1/2097152')


def exact_large_k(orders, c, r):
    u = [Q(1, 4)]*4
    thetas = [u, [Q(17, 64), Q(15, 64), Q(16, 64), Q(16, 64)],
              [Q(15, 64), Q(15, 64), Q(17, 64), Q(17, 64)]]
    tests = 0
    for m in [128, 129, 256]:
        k, j = 8*m, 5*m
        for theta in thetas:
            a = rank_kernel(orders, theta, j, k)
            err = madd(a, scale(-1, c))
            assert max(abs(x) for row in err for x in row) < Q(1, 2**17)
            p = stationary(a)
            assert all(Q(15, 64) < x < Q(17, 64) for x in p)
            assert [x-y for x, y in zip(p, u)] == vecmul(vecmul(p, err), r)
            tests += 1
        # Direct exact endpoint evaluations supplement the uniform bounds.
        assert binomial_tail(k, j, Q(17, 32)) < Q(1, 2**18)
        assert 1-binomial_tail(k, j, Q(47, 64)) < Q(1, 2**18)
        assert k*binomial_atom(k-1, j-1, Q(17, 32)) < Q(1, 32)
        assert k*binomial_atom(k-1, j-1, Q(47, 64)) < Q(1, 32)
    print(f'EXACT LARGE k: {tests} stationary solves at m=128,129,256 verify normalization and perturbation identity')


def main():
    fidelity()
    sampling_crosscheck()
    orders, c, r = cone_and_cycle()
    uniform_arithmetic()
    exact_large_k(orders, c, r)
    print('SCOPE: general attraction uses the stated analytical tail and contraction theorems; no Lean proof performed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
