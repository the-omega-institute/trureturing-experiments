#!/usr/bin/env python3
"""Exact checks for arXiv:2606.18173v1, Algorithm 3.4 / Remark 3.10.

Finite checks are evidence for the proposed general argument, not a proof of
universal reachability. All scalars are integers or Gaussian rationals.
"""
from itertools import combinations, product
import sympy as S

s = S.Symbol("s", real=True)
COUNTS = {"steps": 0, "routes": 0}


def exact(x):
    x = S.simplify(x)
    assert not x.atoms(S.Float), x
    return x


def matrix_exact(A):
    return A.applyfunc(exact)


def psd(A):
    """Hermitian PSD iff every principal minor is nonnegative (exact)."""
    assert A == A.H, A
    n = A.rows
    for r in range(1, n + 1):
        for ix in combinations(range(n), r):
            d = exact(A.extract(ix, ix).det())
            assert d.is_nonnegative is True, (A, ix, d)


def contacts(A, B):
    basis = [v for a in A for v in (B - a).nullspace()]
    return S.Matrix.hstack(*basis).rank() if basis else 0


def solve(P, e):
    u, params = P.gauss_jordan_solve(e)
    u = matrix_exact(u.subs({p: 0 for p in params}))
    assert matrix_exact(P * u - e) == S.zeros(P.rows, 1)
    return u


def step(A, B, e):
    """Literal Steps 3--5, with maximality checked independently by a form."""
    assert e != S.zeros(B.rows, 1)
    assert contacts(A, B) < B.rows
    ee = e * e.H
    coeffs = []
    for a in A:
        P = B - a
        psd(P)
        u = solve(P, e)
        q = exact((u.H * P * u)[0])
        assert q.is_positive is True
        lam = exact(1 / q)
        coeffs.append(lam)
        psd(matrix_exact(P - lam * ee))
        # If s > lam, this exact quadratic form is strictly negative.
        assert exact((u.H * (P - s * ee) * u)[0] - (q - s * q**2)) == 0
    lam = min(coeffs)
    out = matrix_exact(B - lam * ee)
    for a in A:
        psd(out - a)
    assert contacts(A, out) > contacts(A, B)
    COUNTS["steps"] += 1
    return out, coeffs, lam


def fidelity():
    A = [S.Matrix([[3, 1], [1, 3]]), S.Matrix([[5, 1], [1, 1]])]
    B = S.Matrix([[8, 2], [2, 8]])
    B1, cs1, t1 = step(A, B, S.Matrix([1, 1]))
    assert cs1 == [S.Integer(3), S.Rational(5, 2)]
    assert B1 == S.Matrix([[11, -1], [-1, 11]]) / 2
    B2, cs2, t2 = step(A, B1, S.Matrix([-1, 3]))
    assert cs2 == [S.Rational(1, 4), S.Rational(1, 2)]
    assert B2 == S.Matrix([[21, 1], [1, 13]]) / 4
    assert [contacts(A, X) for X in [B, B1, B2]] == [0, 1, 2]
    print("FIDELITY Example 3.9: lambda_i=(3,5/2),(1/4,1/2); B2=[[21,1],[1,13]]/4; contacts=0,1,2")


def steer(A, B, M):
    """Test the proposed direction e=(B-M)v on each exact instance."""
    n = B.rows
    assert contacts(A, M) == n
    D = matrix_exact(B - M)
    psd(D)
    original_rank = D.rank()
    length = 0
    while D != S.zeros(n):
        candidates = [v for a in A for v in (M - a).nullspace() if matrix_exact(D * v) != S.zeros(n, 1)]
        assert candidates
        v = candidates[0]
        e = matrix_exact(D * v)
        c = exact((v.H * D * v)[0])
        assert c.is_positive is True
        t = exact(1 / c)
        old_rank = D.rank()
        B1, coeffs, lam = step(A, B, e)
        assert lam == t, (coeffs, t)
        D1 = matrix_exact(B1 - M)
        psd(D1)
        assert D1.rank() == old_rank - 1
        # The contact vector certifies that any larger update crosses a constraint.
        assert any(matrix_exact((M - a) * v) == S.zeros(n, 1) for a in A)
        assert exact((v.H * (D - s * e * e.H) * v)[0] - (c - s * c**2)) == 0
        B, D = B1, D1
        length += 1
        assert length <= n
    assert B == M
    assert contacts(A, B) == n
    assert length == original_rank
    COUNTS["routes"] += 1
    return length


def scalar_grid():
    cases = 0
    for k in (1, 2, 3):
        for vals in product(range(-2, 3), repeat=k):
            A = [S.Matrix([[a]]) for a in vals]
            M = S.Matrix([[max(vals)]])
            for r in range(4):
                steer(A, M + S.Matrix([[r]]), M)
                cases += 1
    print(f"SCALAR_GRID: n=1; k=1,2,3; all A_i in -2..2; B0=max(A_i)+r, r=0..3; routes={cases}")


def matrix_grid(field):
    cases, singular, complex_cases = 0, 0, 0
    zs = (0, 1, -1) if field == "R" else (S.I, 1 + S.I, -S.I)
    for n in (2, 3):
        diags = list(product((0, 1, 2), repeat=n)) if n == 2 else list(product((0, 1), repeat=n)) + [(2, 2, 2), (1, 2, 1)]
        for z, weighted, mixed, diag in product(zs, (False, True), (False, True), diags):
            T = S.eye(n)
            T[0, 1] = z
            if n == 3:
                T[1, 2], T[0, 2] = 1, -1
            H = []
            for i in range(n):
                h = S.diag(*[0 if j == i else (j + 1 if weighted else 1) for j in range(n)])
                H.append(matrix_exact(T.H * h * T))
            M = S.diag(*[j - 1 for j in range(n)])
            M[0, 1], M[1, 0] = z, S.conjugate(z)
            A = [matrix_exact(M - h) for h in H]
            R = S.eye(n)
            if mixed:
                for j in range(n - 1):
                    R[j + 1, j] = 1 if field == "R" else S.I
            D = matrix_exact(R * S.diag(*diag) * R.H)
            B = matrix_exact(M + D)
            for h in H:
                psd(h)
            for a in A:
                psd(B - a)
            steer(A, B, M)
            cases += 1
            singular += D.rank() < n
            complex_cases += any(exact(S.im(x)) != 0 for x in list(B) + sum((list(a) for a in A), []))
    print(f"MATRIX_GRID {field}: n=k=2,3; all declared contact/residual grid cases; routes={cases}; singular_initial_residual={singular}; nonreal_cases={complex_cases}")


def source_count_correction():
    M = 3 * S.eye(2)
    B = S.diag(4, 3)
    A = [S.diag(3, 2), S.Matrix([[2, -1], [-1, 2]])]
    out, coeffs, lam = step(A, B, S.Matrix([1, 0]))
    assert coeffs == [1, 1] and lam == 1 and out == M
    assert [contacts(A, X) for X in [B, M]] == [0, 2]
    assert (B - M).rank() == 1
    print("ANCILLARY Remark 3.8: B0=diag(4,3), M=3I, A1=diag(3,2), A2=[[2,-1],[-1,2]]; lambda_i=1,1; contact dimension 0->2 in one update")


def main():
    fidelity()
    scalar_grid()
    matrix_grid("R")
    matrix_grid("C")
    source_count_correction()
    print(f"TOTAL: target routes={COUNTS['routes']}; checked maximal updates={COUNTS['steps']}")
    print("Evidence scope: exact finite instances and symbolic quadratic-form identities; universal claim is not proved by this program.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
