#!/usr/bin/env python3
"""Exact diagnostics for the fixed-X higher-eigenvalue defect construction.

Python 3 and SymPy 1.14.0. Every matrix entry, pivot, and spectral threshold
is rational or Gaussian rational. No numerical eigensolver is used.
These checks support a proposed general argument; they do not prove it.
"""
import sympy as s

Q = s.Rational
I = s.I


def clean(M):
    return M.applyfunc(s.expand)


def inertia(M):
    """Exact Hermitian inertia by congruence and Schur complements.

    A nonzero diagonal supplies a 1-by-1 pivot. If every diagonal vanishes,
    a nonzero off-diagonal supplies a 2-by-2 pivot of inertia (1,1,0).
    This also handles singular matrices and thresholds that are eigenvalues.
    """
    M = clean(M)
    assert M == M.H
    original = M.rows
    positive = negative = zero = 0
    while M.rows:
        size = M.rows
        pivot = next((i for i in range(size) if M[i, i] != 0), None)
        if pivot is not None:
            order = [pivot] + [i for i in range(size) if i != pivot]
            M = M.extract(order, order)
            d = M[0, 0]
            assert d.is_Rational
            if d > 0:
                positive += 1
            else:
                negative += 1
            v = M[1:, :1]
            M = clean(M[1:, 1:] - v * v.H / d)
            continue
        pair = next(((i, j) for i in range(size) for j in range(i + 1, size)
                     if M[i, j] != 0), None)
        if pair is None:
            zero += size
            break
        i, j = pair
        order = [i, j] + [k for k in range(size) if k not in (i, j)]
        M = M.extract(order, order)
        z = M[0, 1]
        pivot_inverse = s.Matrix([[0, 1 / s.conjugate(z)], [1 / z, 0]])
        E = M[2:, :2]
        M = clean(M[2:, 2:] - E * pivot_inverse * E.H)
        positive += 1
        negative += 1
    assert positive + negative + zero == original
    return positive, negative, zero


def bracket(M, k, lower, upper):
    """Check lower <= decreasing lambda_k(M) <= upper, with 1-based k."""
    assert 1 <= k <= M.rows and lower <= upper
    p, _, z = inertia(M - lower * s.eye(M.rows))
    assert p + z >= k, (k, lower, (p, z))
    p, _, _ = inertia(M - upper * s.eye(M.rows))
    assert p < k, (k, upper, p)


def block(A, X, B):
    return A.row_join(X).col_join(X.H.row_join(B))


def frobenius_squared(M):
    return s.expand(sum(z * s.conjugate(z) for z in M))


def padded(M, scalar, count):
    return s.diag(M, scalar * s.eye(count))


def source_fidelity():
    # Sababheh Lemma 2.1: triangular minor-axis length is |z|.
    # For beta*I + a*E12 the entire support function is a disk's support.
    a, br, bi, c, d, lam = s.symbols('a br bi c d lam', real=True)
    beta = br + I * bi
    T = s.Matrix([[beta, a], [0, beta]])
    u = c - I * d
    H = clean((u * T + s.conjugate(u) * T.H) / 2)
    expected = (lam - c * br - d * bi)**2 - a**2 * (c**2 + d**2) / 4
    poly = H.charpoly()
    assert s.expand(poly.as_expr().subs(poly.gen, lam) - expected) == 0
    assert s.expand(s.trace(T.H * T) - 2 * (br**2 + bi**2)) == a**2
    # c^2+d^2=1 gives supporting eigenvalues Re(u beta) +/- a/2.
    # In larger dimension the other supporting eigenvalues are Re(u beta).
    for n in range(2, 7):
        X = beta * s.eye(n)
        X[0, 1] = a
        Hn = clean((u * X + s.conjugate(u) * X.H) / 2)
        poly = Hn.charpoly()
        assert s.expand(poly.as_expr().subs(poly.gen, lam) -
            (lam - c * br - d * bi)**(n - 2) * expected) == 0
    # Theorem 3.2 fidelity: Phi_0(a E12)=a, with an exact maximizer.
    for a0 in [Q(1, 2), s.Integer(1), s.Integer(2)]:
        X = s.Matrix([[0, a0], [0, 0]])
        A, B = s.diag(a0, 0), s.diag(0, a0)
        M = block(A, X, B)
        assert inertia(M) == (1, 0, 3)
        assert M.eigenvals() == {2 * a0: 1, 0: 3}
        assert A + B == a0 * s.eye(2)
        assert 2 * a0 - a0 == a0
    print('FIDELITY: Lemma 2.1 disk support identity, n=2..6: exact symbolic equality')
    print('FIDELITY: Theorem 3.2 nilpotent norm-defect maximizers: 3 exact cases')


def fixed_x_family():
    # X=beta I+a E12 is invertible, nonnormal, and has delta_2(X)=a.
    # Symbolic supporting eigenvalues above certify the disk of radius a/2;
    # compression to the first two coordinates attains that full disk.
    beta = 3 + 4 * I
    count = terminal = interior = 0
    for n in range(2, 7):
        for a in [s.Integer(1), s.Integer(2)]:
            X = beta * s.eye(n)
            X[0, 1] = a
            assert s.expand(X.det() - beta**n) == 0
            assert clean(X * X.H - X.H * X) == s.diag(a*a, -a*a, *([0]*(n-2)))
            for j in range(1, n):
                m = n - j
                for q in [2, 4, 8]:
                    if m >= 2:
                        # Rational shift of the source's two-edge construction.
                        t0 = s.Integer(q*q) * (a + 26)
                        error = Q(25, 1) / (2*t0-a) + Q(1, q*q)
                        r = t0 + error
                        A0 = s.diag(r+t0, r-t0, *([r]*(m-2)))
                        B0 = s.diag(r-t0, r+t0, *([r]*(m-2)))
                        t = 10 * (2*r+a+6)
                        A, B = padded(A0, t, j), padded(B0, t, j)
                        M = block(A, X, B)
                        assert inertia(M) == (2*n, 0, 0)
                        bracket(A+B, j+1, 2*r, 2*r)
                        bracket(M, 2*j+1, 2*r+a-error, 2*r+a)
                        assert error > 0 and error < Q(1, q)
                        interior += 1
                    else:
                        # Lift a scalar completion; its low defect is zero.
                        eps = Q(1, q)
                        ell = 10 + eps
                        c2 = a*a
                        t = 2 * (6 + ell + c2/eps + q*q)
                        h = t-5
                        eta = c2/(h-ell)
                        A = B = s.diag(5+eps, *([t]*j))
                        M = block(A, X, B)
                        assert inertia(M) == (2*n, 0, 0)
                        bracket(A+B, j+1, 10+2*eps, 10+2*eps)
                        bracket(M, 2*j+1, ell-eta, 10+2*eps)
                        assert eps+eta < 2*eps
                        terminal += 1
                    count += 1
    print(f'FIXED_X_FAMILY: {count} exact PD completions, n=2..6, every 1<=j<n')
    print(f'INTERIOR: {interior} defects in [a-error,a], error<1/q, q=2,4,8')
    print(f'TERMINAL: {terminal} defects in [-eps-eta,0], eps+eta<2/q')
    print('FAMILY: det(X)=(3+4i)^n; [X,X*]=diag(a^2,-a^2,0,...); delta_2(X)=a>0')


def coupled_lifting():
    count = 0
    for n in range(2, 7):
        for j in range(1, n):
            m = n-j
            eps = Q(1, 4)
            Y = s.zeros(m)
            if m == 1:
                Y[0, 0] = 3+4*I
                A0 = B0 = s.Matrix([[5]])
                ell = 10+eps
                partial = 10+2*eps
            else:
                Y[0, 1] = 1
                A0 = s.diag(1, *([0]*(m-1)))
                B0 = s.diag(0, 1, *([0]*(m-2)))
                ell = 2+eps
                partial = 1+2*eps
            Z = s.Matrix(m, j, lambda a,b: Q((a+2*b)%3-1, 3)+I*Q((a+b)%3-1, 2))
            W = s.Matrix(j, m, lambda a,b: Q((2*a+b)%3-1, 2)+I*Q((a+2*b)%3-1, 3))
            V = s.Matrix(j, j, lambda a,b: (a==b)*2+Q((a+b)%3-1, 2)+I*Q((a+2*b)%3-1, 3))
            X = Y.row_join(Z).col_join(W.row_join(V))
            E = s.zeros(m, j).row_join(Z).col_join(W.H.row_join(s.zeros(m, j)))
            c2 = frobenius_squared(E)
            vbound = 1+sum(abs(s.re(z))+abs(s.im(z)) for z in V)
            t = vbound+ell+2*c2/eps+17
            h = t-vbound
            eta = c2/(h-ell)
            assert h > ell and h > c2/eps
            A = padded(A0+eps*s.eye(m), t, j)
            B = padded(B0+eps*s.eye(m), t, j)
            M = block(A, X, B)
            assert inertia(M) == (2*n, 0, 0)
            assert M[:n, n:] == X
            bracket(M, 2*j+1, ell-eta, ell)
            bracket(A+B, j+1, partial, partial)
            count += 1
    print(f'COUPLED_LIFTING: {count} Gaussian-rational cases with cross-compression blocks')
    print('SPECTRAL_RANKS: 2j+1 in the full completion, j+1 in A+B, checked by exact inertia')


def main():
    # Independent known inertias test the pivot engine, including its 2x2 case.
    assert inertia(s.Matrix([[0, 1+I], [1-I, 0]])) == (1, 1, 0)
    assert inertia(s.diag(3, -2, 0)) == (1, 1, 1)
    source_fidelity()
    fixed_x_family()
    coupled_lifting()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
