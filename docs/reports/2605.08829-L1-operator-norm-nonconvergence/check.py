#!/usr/bin/env python3
"""Exact finite checks for the infinite two-point-block Petz counterexample.

Source: S. Bhattacharya, arXiv:2605.08829v1, equation (2), Theorem 1,
Corollary 3, Corollary 4, and Remarks item 2 after Theorem 2.
Only integers and fractions.Fraction are used; no floating point or external packages.
The infinite-dimensional conclusions require the argument in README.md.
"""
from fractions import Fraction as Q
from itertools import product

I = ((Q(1), Q(0)), (Q(0), Q(1)))
E = ((Q(1, 2), Q(1, 2)), (Q(1, 2), Q(1, 2)))
V = (Q(1), Q(-1))


def matmul(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(2)), Q(0))
                       for j in range(2)) for i in range(2))


def transpose(a):
    return tuple(tuple(a[j][i] for j in range(2)) for i in range(2))


def apply(a, x):
    return tuple(sum((a[i][k] * x[k] for k in range(2)), Q(0)) for i in range(2))


def sub(a, b):
    return tuple(tuple(a[i][j] - b[i][j] for j in range(2)) for i in range(2))


def scale(c, a):
    return tuple(tuple(c * a[i][j] for j in range(2)) for i in range(2))


def power(a, n):
    out = I
    for _ in range(n):
        out = matmul(out, a)
    return out


def channel(epsilon):
    return ((1 - epsilon, epsilon), (epsilon, 1 - epsilon))


def weight(j):
    return Q(1, 2 ** (j + 1))


def norm1(x, w=Q(1)):
    return w * sum(map(abs, x), Q(0)) / 2


def norm1_operator(a):
    # Equal coordinate weights cancel; valid also over complex scalars.
    return max(sum((abs(a[i][j]) for i in range(2)), Q(0)) for j in range(2))


def error(j, n):
    epsilon = Q(1, j + 3)
    lam = (1 - 2 * epsilon) ** 2
    return scale(lam ** n, sub(I, E))


def fidelity():
    # Theorem 1 on a genuine finite tracial algebra, B=1, epsilon=1/4.
    # R=phi*=phi, and T=phi^2; the residual L2 and L1 norms equal (1/4)^n.
    p = channel(Q(1, 4))
    t = matmul(p, p)
    for n in range(33):
        d = sub(power(t, n), E)
        assert d == scale(Q(1, 4) ** n, sub(I, E))
        assert norm1_operator(d) == Q(1, 4) ** n
        assert matmul(transpose(d), d) == scale(Q(1, 4) ** (2 * n), sub(I, E))
    # Corollary 3's L1 contraction on an exact rational grid.
    for x in product((Q(k, 3) for k in range(-6, 7)), repeat=2):
        assert norm1(apply(p, x)) <= norm1(x)
    print('FIDELITY: Theorem 1 finite two-point case, epsilon=1/4, n=0..32;')
    print('  L1 residual norm=(1/4)^n and L2 squared residual=(1/4)^(2n)(I-E).')
    print('  Corollary 3 L1 contraction: 169 exact rational inputs.')


def blocks():
    count = 0
    for j in range(81):
        eps = Q(1, j + 3)
        w = weight(j)
        p = channel(eps)
        t = matmul(p, p)
        lam = (1 - 2 * eps) ** 2
        assert 0 < eps < Q(1, 2) and 0 < lam < 1 and w > 0
        assert all(p[i][k] >= 0 for i in range(2) for k in range(2))
        assert apply(p, (Q(1), Q(1))) == (Q(1), Q(1))
        assert transpose(p) == p  # Hilbert adjoint, since equal weights in this block.
        assert matmul(E, p) == E and matmul(p, E) == E
        assert matmul(E, E) == E
        assert apply(t, V) == tuple(lam * v for v in V)
        for n in range(33):
            d = sub(power(t, n), E)
            assert d == error(j, n)
            assert norm1_operator(d) == lam ** n
            assert norm1(tuple(v / w for v in V), w) == 1
            assert norm1(apply(d, tuple(v / w for v in V)), w) == lam ** n
            count += 1
    for nblocks in range(1, 65):
        assert sum((weight(j) for j in range(nblocks)), Q(0)) + Q(1, 2 ** nblocks) == 1
    print(f'BLOCKS: {count} exact residual identities (j=0..80, n=0..32).')
    print('  Positive stochastic symmetric blocks; B=1; Petz R=phi; fixed projection E.')
    print('  Geometric weights: each finite sum plus its exact tail is 1.')


def truncations():
    for n in (0, 1, 2, 4, 8, 16, 32):
        for nblocks in (1, 2, 4, 8, 16, 32, 64, 128):
            computed = max(norm1_operator(error(j, n)) for j in range(nblocks))
            formula = (1 - Q(2, nblocks + 2)) ** (2 * n)
            assert computed == formula
    # A fixed N does converge. For each n we instead take N(n)=8n+1,
    # so the last included block has j=8n and its norm exceeds 1/2.
    for n in range(1, 129):
        nblocks = 8 * n + 1
        value = norm1_operator(error(8 * n, n))
        bernoulli_lower = 1 - Q(4 * n, 8 * n + 3)
        assert value >= bernoulli_lower > Q(1, 2)
        if n in (1, 2, 4, 8, 16):
            print(f'MOVING_TRUNCATION n={n} N={nblocks} norm={value} lower={bernoulli_lower}>1/2')
    # For arbitrary rational eta=1/d and n, this explicit block has error >=1-eta.
    for n in range(1, 65):
        for d in range(2, 17):
            j = 4 * n * d
            value = norm1_operator(error(j, n))
            assert value >= 1 - Q(4 * n, j + 3) > 1 - Q(1, d)
    print('TRUNCATIONS: exact formula ||T_N^n-E_N||_1=(1-2/(N+2))^(2n).')
    print('  n=0,1,2,4,8,16,32; N=1,2,4,8,16,32,64,128.')
    print('  Moving N(n)=8n+1 gives norms >1/2 for n=1..128.')
    print('  Explicit near-1 witnesses: n=1..64, eta=1/d, d=2..16.')
    print('  Infinite norm=1 and strong L1 convergence require the stated tail/limit argument;')
    print('  finite truncations alone do not certify an infinite supremum or a limit.')


def main():
    fidelity()
    blocks()
    truncations()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
