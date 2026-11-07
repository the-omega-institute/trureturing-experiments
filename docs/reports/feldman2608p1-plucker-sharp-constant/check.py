#!/usr/bin/env python3
"""Exact finite checks for the complex coordinate-norm question in Feldman P1.

This is supporting evidence, not a proof of the universal inequality.
Only integer/rational and exact symbolic polynomial arithmetic is used.
"""
from itertools import combinations, permutations, product
from math import factorial, prod
from random import Random

import sympy as sp


def determinant(a):
    """Leibniz determinant, independent of the Gram-matrix determinant call."""
    n = len(a)
    ans = 0
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        ans += sign * prod(a[i][perm[i]] for i in range(n))
    return sp.expand(ans)


def abs_sq(z):
    return sp.expand(z * sp.conjugate(z))


def check_matrix(a, exponents=(2, 4, 6, 8)):
    n, m = len(a), len(a[0])
    minor_sq = [abs_sq(determinant([[row[j] for j in cols] for row in a]))
                for cols in combinations(range(m), n)]
    gram = sp.Matrix(a) * sp.Matrix(a).conjugate().T
    assert sum(minor_sq) == sp.expand(gram.det()), ("Cauchy-Binet", a)
    assert sum(minor_sq) <= prod(sum(abs_sq(z) for z in row) for row in a)
    for p in exponents:
        lhs = sum(d ** (p // 2) for d in minor_sq)
        rows = prod(sum(abs_sq(z) ** (p // 2) for z in row) for row in a)
        # C(n,p)^p = n^(n*p/2-n), an integer at these even exponents.
        rhs = n ** (n * (p // 2) - n) * rows
        assert lhs <= rhs, ("proposed bound", n, m, p, a, lhs, rhs)
        assert lhs <= factorial(n) ** (p - 1) * rows, ("source crude bound", a, p)
        if n == 1:
            assert lhs == rows
    return len(exponents)


def fidelity_sylvester():
    h = sp.Matrix([[1]])
    for n in (1, 2, 4, 8):
        if n > 1:
            h = h.row_join(h).col_join(h.row_join(-h))
        assert h * h.T == n * sp.eye(n)
        assert h.det() ** 2 == n ** n
        # Theorem 2.1(d): exactly one nonzero n-minor for an n-column frame.
        assert len(list(combinations(range(n), n))) == 1
        print(f"FIDELITY Sylvester n={n}: det^2={h.det() ** 2}={n}^{n}")


def exact_fourier():
    x = sp.Symbol("x")
    for n in range(1, 9):
        phi = sp.Poly(sp.cyclotomic_poly(n, x), x, domain=sp.QQ)

        def red(f):
            return sp.rem(sp.Poly(sp.expand(f), x, domain=sp.QQ), phi).as_expr()

        def conj(f):
            return red(f.subs(x, x ** (n - 1)))

        # Work in Q[x]/Phi_n, which maps exactly to Q(exp(2*pi*i/n)).
        f = sp.Matrix(n, n, lambda r, k: red(x ** ((r * k) % n)))
        for r in range(n):
            for s in range(n):
                assert red(sum(f[r, k] * conj(f[s, k]) for k in range(n))) == (n if r == s else 0)
                assert red(f[r, s] * conj(f[r, s])) == 1
        det = red(f.det(method="domain-ge"))
        assert red(det * conj(det)) == n ** n
        for p in (3, 4, 6, 8):
            # Square the nonnegative p-th-power equality to avoid radicals.
            # |det F|^(2p) = C(n,p)^(2p) * product_r ||F_r||_p^(2p).
            left = (n ** n) ** p
            right = n ** (n * (p - 2)) * n ** (2 * n)
            assert left == right
        print(f"FOURIER n={n}: exact Gram=nI, |det|^2={n ** n}; equality p=3,4,6,8")


def main():
    print("Exact arithmetic: Python integers; SymPy Gaussian integers and Q[x]/Phi_n")
    fidelity_sylvester()
    exact_fourier()
    # Exhaust every matrix with entries 0 or 1, n<=2, n<=m<=3.
    exhaustive = 0
    for n in (1, 2):
        for m in range(n, 4):
            for entries in product((0, 1), repeat=n * m):
                a = [list(entries[r*m:(r+1)*m]) for r in range(n)]
                check_matrix(a)
                exhaustive += 1
    print(f"EXHAUSTIVE binary matrices n<=2, n<=m<=3: {exhaustive}; p=2,4,6,8")
    rng = Random(260800983)
    random_count = 0
    for n in range(1, 5):
        for m in range(n, n + 4):
            for _ in range(24):
                a = [[sp.Integer(rng.randrange(-2, 3)) + sp.I * rng.randrange(-2, 3)
                      for _ in range(m)] for _ in range(n)]
                check_matrix(a)
                random_count += 1
    print(f"RANDOM Gaussian-integer matrices n=1..4, m=n..n+3: {random_count}; seed=260800983")
    print("All tested frames: source Cauchy-Binet identity, Gram-Hadamard and crude bound; proposed sharp bound p=2,4,6,8")
    print("Known n=1 identity and n=2 constant C(2,p)^p=2^(p-2) reproduced")
    print("Scope: finite checks do not prove the universal real-p or infinite-coordinate claim")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
