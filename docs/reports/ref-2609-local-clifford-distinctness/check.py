#!/usr/bin/env python3
"""Exact source-fidelity and settlement checks for arXiv:2609.03065v1 I.C.a.

Evidence uses integers, fractions.Fraction and exact SymPy algebra only.
The collision check reproduces the source's Proposition III.5. This file
contains no Lean proof and does not check cryptographic security.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import sympy as s

I2 = s.eye(2)
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -s.I], [s.I, 0]])
Z = s.diag(1, -1)
H = (X + Z) / s.sqrt(2)
S = s.diag(1, s.I)


def clean(A):
    return A.applyfunc(s.simplify)


def canon(U):
    pivot = next(x for x in U if x != 0)
    return tuple(s.simplify(x / pivot) for x in U)


def reps():
    six = [I2, H, S, H*S, S*H, H*S*H]
    result = [clean(X**a * Z**b * V)
              for a, b, V in product(range(2), range(2), six)]
    keys = {canon(U) for U in result}
    assert len(keys) == 24
    # H and S generate Clifford; closure and 24 distinct classes check the list.
    for U in result:
        assert clean(U.conjugate().T * U) == I2
        assert canon(H*U) in keys and canon(S*U) in keys
    return result


def epsilon(n, t):
    N = 2**n
    assert n >= 2 and 0 <= t <= N//2
    return 4*(comb(t, 2)*F(2, 3)**n + F(t, N)) + F(t*(t-1), N-t+1)


def frac(x):
    return F(int(s.numer(x)), int(s.denom(x)))


def dnp_t2(n):
    """All eigenvalues of the averaged DNP projector, without floating point.

    On Dist, DNP consists of zero-row-sum and zero-column-sum matrices.
    Remove the last column equation to obtain independent constraints A.
    Its projector is I-A^T(AA^T)^(-1)A. The two-copy local Clifford twirl
    is the Hilbert-Schmidt projection onto span{I,SWAP}; the inverse Gram
    matrix is [[1/3,-1/6],[-1/6,1/3]]. Commuting local swaps give the full
    spectrum, hence a universal check over all density matrices at t=2.
    """
    N = 2**n
    edges = [(x, y) for x in range(N) for y in range(N) if x != y]
    A = s.zeros(2*N-1, len(edges))
    for j, (x, y) in enumerate(edges):
        A[x, j] = 1
        if y < N-1:
            A[N+y, j] = 1
    Ginv = (A*A.T).inv()
    assert (A*A.T)*Ginv == s.eye(2*N-1)
    Ginv = [[frac(v) for v in row] for row in Ginv.tolist()]

    def support(x, y):
        return [x] + ([N+y] if y < N-1 else [])

    def entry(x, y, u, v):
        if x == y or u == v:
            return F(0)
        return F(int(x == u and y == v)) - sum(
            (Ginv[i][j] for i in support(x, y) for j in support(u, v)), F(0))

    traces = []
    for mask in range(2**n):
        value = F(0)
        for x, y in edges:
            delta = (x ^ y) & mask
            value += entry(x ^ delta, y ^ delta, x, y)
        traces.append(value)
    assert traces[0] == N*(N-1) - (2*N-1)
    coeffs = [sum((traces[b] * F((-1)**((a ^ b).bit_count()),
                                3**n * 2**((a ^ b).bit_count()))
                   for b in range(2**n)), F(0)) for a in range(2**n)]
    eigenvalues = [sum((coeffs[a] * (-1)**((a & negative).bit_count())
                       for a in range(2**n)), F(0))
                   for negative in range(2**n)]
    assert min(eigenvalues) >= 0 and max(eigenvalues) <= 1
    assert min(eigenvalues) >= 1-epsilon(n, 2)
    print(f'DNP universal t=2: n={n}, min_overlap={min(eigenvalues)}, '
          f'epsilon={epsilon(n, 2)}, nonvacuous={epsilon(n, 2) < 1}')


def main():
    C = reps()
    print('Clifford classes modulo scalar phase: 24; unitarity and H/S closure exact')
    counts = {(j, sign): 0 for j in range(3) for sign in [-1, 1]}
    for U in C:
        image = clean(U.conjugate().T * Z * U)
        hit = [(j, sign) for j, P in enumerate([X, Y, Z]) for sign in [-1, 1]
               if image == sign*P]
        assert len(hit) == 1
        counts[hit[0]] += 1
    assert set(counts.values()) == {4}
    print('Z conjugation: each of +/-X, +/-Y, +/-Z has multiplicity 4')
    pairs = list(product(range(2), repeat=2))
    swap = s.Matrix([[int(a == d and b == c) for c, d in pairs] for a, b in pairs])
    Qeq = s.diag(1, 0, 0, 1)
    twirl = clean(sum((s.kronecker_product(U, U).conjugate().T * Qeq *
                       s.kronecker_product(U, U) for U in C), s.zeros(4))/24)
    assert twirl == (s.eye(4)+swap)/3
    assert twirl.eigenvals() == {s.Rational(2, 3): 3, s.Integer(0): 1}
    print('Fidelity, Proposition III.5: twirl=(I+SWAP)/3; spectrum 2/3 (3), 0 (1)')
    plus = s.ones(2)/2
    first = clean(sum((U.conjugate().T * plus * U for U in C), s.zeros(2))/24)
    assert first == I2/2
    print('Uniform-projector first moment: I/2 exactly')
    for a, b in product(range(4), repeat=2):
        E = s.zeros(4)
        E[a, b] = 1
        actual = clean(sum((s.kronecker_product(U, U).conjugate().T * E *
                            s.kronecker_product(U, U) for U in C), s.zeros(4))/24)
        expected = ((s.trace(E)/3 - s.trace(swap*E)/6)*s.eye(4) +
                    (s.trace(swap*E)/3 - s.trace(E)/6)*swap)
        assert actual == expected
    print('Full two-copy Clifford twirl: all 16 matrix units checked exactly')
    for n in range(1, 5):
        T = s.ImmutableSparseMatrix(twirl)
        for _ in range(n-1):
            T = s.ImmutableSparseMatrix(s.kronecker_product(T, twirl))
        lam = s.Rational(2, 3)**n
        assert T.T == T and T*T == lam*T
        e0 = s.zeros(4**n, 1)
        e0[0] = 1
        assert T*e0 == lam*e0
        assert max(sum(abs(x) for x in T.row(i)) for i in range(T.rows)) == lam
        print(f'Tensor collision: n={n}, dimension={4**n}, exact norm={lam}')
    for n in range(2, 5):
        N = 2**n
        for t in range(2, min(4, N//2)+1):
            checked = 0
            for xs in product(range(N), repeat=t):
                collisions = sum(xs[i] == xs[j] for i in range(t) for j in range(i+1, t))
                assert int(len(set(xs)) < t) <= collisions
                checked += 1
            print(f'Commuting collision union inequality: n={n}, t={t}, tuples={checked}')
    for n in range(2, 6):
        dnp_t2(n)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
