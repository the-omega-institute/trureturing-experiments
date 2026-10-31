#!/usr/bin/env python3
"""Exact checks for Arvind--Hunter-Jones, arXiv:2610.06754v1, Conjecture 1.

Requires Python 3 and SymPy. No floating point numbers are used as evidence.
The unbounded replica claim is supported by the general argument in the
preregistration; the finite replica checks here are independent diagnostics.
"""

from itertools import permutations, product
import sympy as S

I = S.I
Id = S.eye(2)
X = S.Matrix([[0, 1], [1, 0]])
Y = S.Matrix([[0, -I], [I, 0]])
Z = S.diag(1, -1)
paulis = [Id, X, Y, Z]


def tensor(*args):
    out = S.ones(1)
    for arg in args:
        out = S.kronecker_product(out, arg)
    return out


def equal(A, B):
    assert A.shape == B.shape
    assert all(S.cancel(a - b) == 0 for a, b in zip(A, B))


def bracket(A, B):
    return (A * B - B * A) / (2 * I)


def core(t):
    return S.Matrix([[1/t, 0, 0, 0], [0, 0, -I*t, 0],
                     [0, -I*t, 0, 0], [0, 0, 0, 1/t]])


def reshuffle(U):
    # Row (a,c), column (b,d); U has output (c,d), input (a,b).
    return S.Matrix(4, 4, lambda r, s:
                    U[2*(r % 2) + s % 2, 2*(r // 2) + s // 2])


def partial_trace_right(A):
    return S.Matrix(2, 2, lambda r, s: sum(A[2*r+b, 2*s+b] for b in range(2)))


def symbolic_checks():
    t = S.Symbol('t', nonzero=True)
    U = core(t)
    # Unit-circle involution: conjugate(t)=1/t.
    Udag = S.Matrix([[t, 0, 0, 0], [0, 0, I/t, 0],
                    [0, I/t, 0, 0], [0, 0, 0, t]])
    equal(U * Udag, S.eye(4))
    equal(reshuffle(U) * reshuffle(Udag).T, S.eye(4))
    # Source fidelity: dual unitarity spreads every traceless one-site operator.
    for P in [X, Y, Z]:
        equal(partial_trace_right(U * tensor(P, Id) * Udag), S.zeros(2))
    print('FIDELITY: symbolic dual unitarity and traceless one-site spreading: OK')

    s = (t**2 - t**-2)/(2*I)
    c = (t**2 + t**-2)/2
    A = s*tensor(X, Id) + c*tensor(Y, Z)
    equal(Udag * tensor(Id, X) * U, A)
    equal(Udag * tensor(Id, Y) * U, s*tensor(Y, Id) - c*tensor(X, Z))
    equal(Udag * tensor(Id, Z) * U, tensor(Z, Id))
    equal(bracket(tensor(Id, Y), A), c*tensor(Y, X))
    equal(bracket(tensor(Id, Z), tensor(Y, X)), tensor(Y, Y))
    equal(bracket(tensor(Id, X), tensor(Y, Y)), tensor(Y, Z))
    equal((A-c*tensor(Y, Z))/s, tensor(X, Id))
    equal(bracket(tensor(Z, Id), tensor(X, Id)), tensor(Y, Id))
    for P in [X, Y, Z]:
        equal(bracket(tensor(Z, Id), tensor(Y, P)), -tensor(X, P))
        equal(bracket(tensor(X, Id), tensor(Y, P)), tensor(Z, P))
    basis = [tensor(P, Q) for j, P in enumerate(paulis)
             for k, Q in enumerate(paulis) if (j, k) != (0, 0)]
    gram = S.Matrix(15, 15, lambda j, k: S.trace(basis[j].H * basis[k]))
    equal(gram, 4*S.eye(15))
    print('SETTLEMENT: all parameter-dependent conjugations, generating brackets, and Pauli Gram: OK')


def perm_matrix(m, pi):
    configs = list(product(range(2), repeat=m))
    index = {x: j for j, x in enumerate(configs)}
    inv = [pi.index(j) for j in range(m)]
    P = S.zeros(2**m)
    for j, x in enumerate(configs):
        P[index[tuple(x[inv[k]] for k in range(m))], j] = 1
    return P


def replica(U, m):
    configs = list(product(range(2), repeat=m))
    pairs = list(product(configs, configs))
    return S.Matrix(len(pairs), len(pairs), lambda r, s:
                    S.prod(U[2*pairs[r][0][k]+pairs[r][1][k],
                             2*pairs[s][0][k]+pairs[s][1][k]] for k in range(m)))


def vec(A):
    return S.Matrix(list(A))


def unvec(v, n):
    return S.Matrix(n, n, list(v))


def sector_basis(m):
    n = 2**m
    perms = [perm_matrix(m, pi) for pi in permutations(range(m))]
    diag = S.Matrix.hstack(*[vec(tensor(P, P)) for P in perms])
    diag = S.Matrix.hstack(*diag.columnspace())
    gram_inv = (diag.H*diag).inv()
    cols = []
    for a, b in product(range(n), repeat=2):
        E = S.zeros(n)
        E[a, b] = 1
        for P in perms:
            v = vec(tensor(E, P))
            cols.append(v - diag*gram_inv*(diag.H*v))
    B = S.Matrix.hstack(*cols)
    B = S.Matrix.hstack(*B.columnspace())
    equal(diag.H*B, S.zeros(diag.cols, B.cols))
    assert B.cols == {1: 3, 2: 30}[m]
    return perms, diag, B


def positive_definite_exact(A):
    equal(A, A.H)
    # Exact Hermitian Schur complements; positive pivots certify positive definiteness.
    pivots = []
    while A.rows:
        pivot = S.cancel(A[0, 0])
        assert pivot.is_Rational and pivot > 0, pivot
        pivots.append(pivot)
        if A.rows == 1:
            break
        A = (A[1:, 1:] - A[1:, :1]*A[:1, 1:]/pivot).applyfunc(S.cancel)
    return min(pivots)


def small_replica_checks():
    q = S.Rational
    phases = [(q(3, 5)+I*q(4, 5)), (q(5, 13)+I*q(12, 13))]
    # A genuinely complex exact SU(2) dressing, used on all four legs.
    D = S.Matrix([[q(1, 2)+I*q(1, 2), q(1, 2)+I*q(1, 2)],
                  [-q(1, 2)+I*q(1, 2), q(1, 2)-I*q(1, 2)]])
    equal(D.H*D, S.eye(2))
    assert D.det() == 1
    gates = [(f'phase-{j}-bare', core(t)) for j, t in enumerate(phases)]
    gates.append(('phase-0-local-dressing', tensor(D, D.H)*core(phases[0])*tensor(D.H, D)))
    gates = [(label, U.applyfunc(S.expand_complex)) for label, U in gates]
    for m in [1, 2]:
        perms, diag, B = sector_basis(m)
        Q = B.H*B
        Qinv = Q.inv()
        for label, U in gates:
            equal(U.H*U, S.eye(4))
            R = reshuffle(U)
            equal(R.H*R, S.eye(4))
            Um = replica(U, m).applyfunc(S.expand_complex)
            for P in perms:
                equal(Um*tensor(P, P)*Um.H, tensor(P, P))
            images = S.Matrix.hstack(*[vec(Um*unvec(B[:, j], 2**(2*m))*Um.H)
                                      for j in range(B.cols)])
            T = (B.H*images).applyfunc(S.cancel)
            if m == 1:
                equal(T, S.zeros(3))
                print(f'FIDELITY: {label}, m=1: compression zero, source gap 1: OK')
            else:
                H = (Q - T.H*Qinv*T).applyfunc(S.cancel)
                pivot = positive_definite_exact(H)
                print(f'SETTLEMENT: {label}, m=2, dim(B)=30: strict contraction; '
                      f'min exact Schur pivot={pivot}: OK')
    # Dependence among single-qubit permutation operators must not be ignored.
    P3 = S.Matrix.hstack(*[vec(perm_matrix(3, pi)) for pi in permutations(range(3))])
    assert P3.rank() == 5
    print('FIDELITY: m=3: six permutation operators have span dimension 5: OK')


if __name__ == '__main__':
    symbolic_checks()
    small_replica_checks()
    print('ALL_CHECKS_OK')
