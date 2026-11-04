#!/usr/bin/env python3
"""Exact checks for arXiv:2610.34388v1, Section 11.2(5).

Uses only integers and exact rational complex algebra. Finite checks support,
but do not prove, the proposed all-dimensional theorem.
"""
from itertools import product
import sympy as s

I, Q = s.I, s.Rational


def clean(x):
    return s.expand_complex(s.cancel(s.expand(x)))


def cm(M):
    return M.applyfunc(clean)


def eq(M, N):
    assert cm(M - N) == s.zeros(*M.shape), (M, N)


def realpart(X):
    return cm((X + X.H) / 2)


def positive(x):
    x = clean(x)
    assert x.is_Rational and x > 0, x


def pd(M):
    """Exact Sylvester criterion for the Hermitian matrix passed in."""
    eq(M, M.H)
    for k in range(1, M.rows + 1):
        positive(M[:k, :k].det())


def fidelity():
    # Proposition 5.44(d): reproduce every displayed exact quantity.
    P = s.Matrix([[2, 19, 20], [19, 637, 491], [20, 491, 416]])
    U = s.Matrix([[9965, -1180, 3040], [2020, 9907, -2776],
                  [-2560, 3224, 9643]]) / 10485
    rho = s.diag(10000, 1, 1)
    eq(U.H * U, s.eye(3))
    pd(P)
    X = U * P
    eq(X.H * X, P * P)  # P positive definite, hence P = |X|.
    A = realpart(X)
    pd(A)
    pminors = [P[:k, :k].det() for k in range(1, 4)]
    aminors = [clean(A[:k, :k].det()) for k in range(1, 4)]
    assert pminors == [2, 913, 16006]
    assert aminors == [Q(11662, 2097), Q(6082302974, 109935225),
                       Q(36535906, 7329015)]
    num, den = clean(s.trace(rho * A)), s.trace(rho * P)
    assert num == Q(39575293, 699) and den == 21053
    ratio = clean(num / den)
    assert ratio == Q(39575293, 14716047) and ratio > Q(268, 100)
    print('FIDELITY Proposition 5.44(d): P minors =', pminors)
    print('FIDELITY Re(X) minors =', aminors)
    print('FIDELITY quotient =', ratio, '> 67/25')


def family(n, t):
    # A rational-coordinate version of the square-root triangular family.
    # T_{jk} = 2*c_k above the diagonal; W = diag(1/c_j).
    # (T*W + W*T.H)/2 = ee.H. All entries lie in Q(i).
    c = [clean((1 - j*j*t*t)/(1 + j*j*t*t)) for j in range(1, n+1)]
    lam = [clean((1 + I*j*t)/(1 - I*j*t)) for j in range(1, n+1)]
    assert 0 < t < Q(1, n) and len(set(lam)) == n
    for x in c:
        positive(x)
    D, W = s.diag(*lam), s.diag(*[1/x for x in c])
    T = s.zeros(n)
    for j in range(n):
        T[j, j] = lam[j]
        for k in range(j+1, n):
            T[j, k] = 2*c[k]
    e = s.ones(n, 1)
    eq((T*W + W*T.H)/2, e*e.H)
    eq(D.H*D, s.eye(n))
    S = s.eye(n)
    for j in range(n):
        for i in reversed(range(j)):
            S[i, j] = clean(-sum(T[i, k]*S[k, j]
                                for k in range(i+1, j+1))/(lam[i]-lam[j]))
    eq(T*S, S*D)
    R = cm(S.inv())
    P = cm(R*W*R.H)
    X = cm(D*P)
    A = realpart(X)
    y = cm(R*e)
    eq(A, y*y.H)  # Exact positive semidefinite rank-one certificate.
    assert A.rank() == 1
    pd(P)
    assert X.det() != 0
    eq(X.H*X, P*P)  # Thus |X| = P by uniqueness of the positive square root.
    pd(n*P-A)
    L, C = sum(c), n-Q(1, 100)
    assert n > L > C > 0
    z = cm(R.H.inv()*s.Matrix(c))
    eq(z.H*A*z, s.Matrix([[L*L]]))
    eq(z.H*P*z, s.Matrix([[L]]))
    b = clean(s.trace(P))
    eta = clean(L*(L-C)/(2*C*b))
    positive(eta)
    rho = cm(z*z.H + eta*s.eye(n))
    pd(rho)
    num, den = clean(s.trace(rho*A)), clean(s.trace(rho*P))
    positive(den)
    ratio = clean(num/den)
    assert ratio > C and ratio < n
    assert ratio >= 2*C*L/(C+L) > C
    B = cm((X-X.H)/(2*I))
    eq(B, B.H)
    eq(A+I*B, X)
    # A singular accretive extension checks the support/rank boundary.
    X0, P0, A0 = s.diag(X, 0), s.diag(P, 0), s.diag(A, 0)
    assert X0.rank() == n and X0[:, n] == s.zeros(n+1, 1)
    eq(X0.H*X0, P0*P0)
    eq(realpart(X0), A0)
    print(f'FAMILY n={n}, t={t}: exact identities, positivity, rank support; '
          f'{C} < quotient < {n}')


def small_grid():
    # Exhaust the stated 2x2 grid. Positivity and the order bound are checked
    # by all principal minors, without numerical eigenvalues.
    units = [s.Integer(1), I, Q(3, 5)+I*Q(4, 5), Q(5, 13)+I*Q(12, 13)]
    tested = accepted = 0
    for p, q, a, b in product(range(1, 3), range(1, 3), range(-1, 2), range(-1, 2)):
        P = s.Matrix([[p, a+I*b], [a-I*b, q]])
        if P.det() <= 0:
            continue
        for u, v in product(units, repeat=2):
            tested += 1
            X = s.diag(u, v)*P
            A = realpart(X)
            if A[0, 0] < 0 or A[1, 1] < 0 or clean(A.det()) < 0:
                continue
            accepted += 1
            H = cm(2*P-A)
            assert H[0, 0] >= 0 and H[1, 1] >= 0 and clean(H.det()) >= 0
            eq(X.H*X, P*P)
    assert accepted > 0
    print(f'GRID 2x2: {tested} inputs with P>0; {accepted} accretive; '
          'all satisfy Re(X) <= 2|X| exactly')


def main():
    print('ARITHMETIC: exact integers and SymPy rationals over Q(i); SymPy', s.__version__)
    fidelity()
    for n in range(1, 6):
        family(n, Q(1, 100*n*n))
    small_grid()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
