#!/usr/bin/env python3
"""Exact checks for BCEHK, arXiv:2206.15410v1, Conjectures 14 and 20.

All numerical certificates use integers or SymPy rationals and rational real-root
isolation. U has columns e_i-e_(n-1); (U.T*U)^-1 U.T*S*U represents the orthogonal
compression of S to the zero-sum space in this nonorthonormal basis.
"""
import hashlib
from itertools import product
from pathlib import Path
import sympy as sp

R = sp.Rational


def laplacian(A):
    n = A.rows
    assert A.cols == n
    assert all(A[i, i] == 0 for i in range(n))
    assert all(a in (0, 1) for a in A)
    return sp.diag(*[sum(A.row(i)) for i in range(n)]) - A


def compression(L):
    n = L.rows
    U = sp.Matrix.vstack(sp.eye(n - 1), -sp.ones(1, n - 1))
    G = U.T * U
    S = (L + L.T) / 2
    return G.inv() * U.T * S * U


def root_intervals(L):
    p = compression(L).charpoly().as_poly()
    intervals = p.intervals(eps=R(1, 10**8))
    assert sum(multiplicity for _, multiplicity in intervals) == L.rows - 1
    return intervals


def rayleigh(L, v):
    v = sp.Matrix(v)
    assert sum(v) == 0 and (v.T * v)[0] > 0
    return (v.T * L * v)[0] / (v.T * v)[0]


def fidelity():
    # Corollary 3 in the rendered v1 (cor:cycle_spread in TeX).
    for n, expected in [(2, 0), (3, 0), (4, 1), (6, R(3, 2))]:
        A = sp.zeros(n)
        for i in range(n):
            A[i, (i + 1) % n] = 1
        vals = compression(laplacian(A)).eigenvals()
        assert all(v.is_Rational for v in vals)
        actual = max(vals) - min(vals)
        assert actual == expected
        print(f'FIDELITY directed cycle n={n}: spread={actual}')


def balanced_small():
    expected_counts = {2: 2, 3: 10, 4: 152}
    for n in (2, 3, 4):
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        count = zero = positive = 0
        min_positive_lower_bound = None
        for bits in product((0, 1), repeat=len(pairs)):
            outgoing = [0] * n
            incoming = [0] * n
            for (i, j), bit in zip(pairs, bits):
                outgoing[i] += bit
                incoming[j] += bit
            if outgoing != incoming:
                continue
            A = sp.zeros(n)
            for (i, j), bit in zip(pairs, bits):
                A[i, j] = bit
            L = laplacian(A)
            S = (L + L.T) / 2
            assert S * sp.ones(n, 1) == sp.zeros(n, 1)
            # Check the degree Rayleigh identity independently on every graph.
            for i in range(n):
                p = -sp.ones(n, 1) / n
                p[i] += 1
                assert rayleigh(L, p) == R(n * outgoing[i], n - 1)
            intervals = root_intervals(L)
            (alo, ahi), _ = intervals[0]
            (blo, bhi), _ = intervals[-1]
            if len(intervals) == 1:
                assert alo == ahi
                zero += 1
            else:
                lower = blo - ahi
                assert lower >= 1, (n, A, intervals)
                positive += 1
                if min_positive_lower_bound is None or lower < min_positive_lower_bound:
                    min_positive_lower_bound = lower
            count += 1
        assert count == expected_counts[n]
        print(f'C14 n={n}: labeled_balanced={count}, zero={zero}, '
              f'positive={positive}, positive_spread_lower_bound={min_positive_lower_bound}')
    # Check the algebraic variance identity and strict estimate in a larger
    # finite parameter window; this is not a proof of the universal inequality.
    checked = 0
    for n in range(5, 65):
        for r in range(2, n - 2, 2):
            mu = -R(r, n - 1)
            variance = R(n * r - r * r, n - 1) - mu * mu
            assert variance == R(n * r * (n - 1 - r), (n - 1)**2)
            assert variance >= R(2 * n * (n - 3), (n - 1)**2) > 1
            checked += 1
    print(f'C14 variance parameters: {checked} exact checks, 5<=n<=64')


def slope_family():
    # Figure 3, ordered left, upper middle, lower middle, upper right, lower right.
    A = sp.Matrix([[0, 0, 1, 0, 0],
                   [1, 0, 1, 1, 1],
                   [1, 1, 0, 1, 1],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0]])
    L = laplacian(A)
    z = sp.Symbol('z')
    p = compression(L).charpoly(z).as_poly()
    assert sp.expand(p.as_expr() - z * (20*z**3 - 180*z**2 + 437*z - 190)/20) == 0
    intervals = p.intervals(eps=R(1, 10**10))
    (alo, ahi), _ = intervals[0]
    (blo, bhi), _ = intervals[-1]
    assert alo == ahi == 0
    # Reproduce the source's rounded spread 5.035 by a rational enclosure.
    assert R(10069, 2000) < blo <= bhi < R(10071, 2000)
    print(f'FIDELITY Figure 3: alpha=0, beta in [{blo}, {bhi}], rounds to 5.035')
    print(f'C20 base adjacency: {A.tolist()}')
    print(f'C20 compressed characteristic polynomial: {p.as_expr()}')
    x = [2, 26, -26, -1, -1]
    y = [0, 0, 0, 1, -1]
    qx, qy = rayleigh(L, x), rayleigh(L, y)
    assert qx == R(3408, 679) and qy == 0
    assert qx - qy - 5 == R(13, 679) > 0
    print(f'C20 x={x}, y={y}, Rayleigh(x)={qx}, Rayleigh(y)={qy}, excess={qx-qy-5}')
    for m in range(1, 9):
        Am = sp.kronecker_product(A, sp.ones(m))
        Lm = laplacian(Am)
        assert Lm == m * sp.kronecker_product(sp.diag(1, 4, 4, 0, 0), sp.eye(m)) - Am
        xm = [value for value in x for _ in range(m)]
        ym = [value for value in y for _ in range(m)]
        assert rayleigh(Lm, xm) == m * qx
        assert rayleigh(Lm, ym) == 0
        excess = rayleigh(Lm, xm) - rayleigh(Lm, ym) - 5*m
        assert excess == R(13*m, 679) > 0
        print(f'C20 independent-cluster blowup m={m}: order={5*m}, '
              f'spread>={m*qx}, spread-order>={excess}')


if __name__ == '__main__':
    print(f'Python exact arithmetic; SymPy {sp.__version__}')
    print('check.py sha256=' + hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    fidelity()
    balanced_small()
    slope_family()
    print('ALL_CHECKS_OK')
