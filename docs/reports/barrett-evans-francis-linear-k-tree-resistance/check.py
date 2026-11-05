#!/usr/bin/env python3
"""Exact checks for the straight linear k-tree resistance increment.

Python 3 standard library only.  All evidence uses integers or Fraction.
Grounded Kirchhoff equations and interval-Gram equations are assembled
independently.  These finite checks do not establish an infinite limit.
"""
from fractions import Fraction as F


def solve(matrix, rhs):
    """Rational Gaussian elimination with row pivoting; no approximate arithmetic."""
    n = len(rhs)
    a = [[F(v) for v in row] + [F(rhs[i])] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [v / scale for v in a[col]]
        for row in range(col + 1, n):
            scale = a[row][col]
            if scale:
                a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
    x = [F(0)] * n
    for i in reversed(range(n)):
        x[i] = a[i][n] - sum((a[i][j] * x[j] for j in range(i + 1, n)), F(0))
    return x


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), F(0))


def matvec(a, x):
    return [dot(row, x) for row in a]


def laplacian(n, k):
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if 0 < abs(i - j) <= k:
                a[i][i] += 1
                a[j][j] += 1
                a[i][j] -= 1
                a[j][i] -= 1
    return a


def grounded_resistance(n, k):
    a = laplacian(n, k)
    v = solve([row[:-1] for row in a[:-1]], [1] + [0] * (n - 2)) + [F(0)]
    assert matvec(a, v) == [F(1)] + [F(0)] * (n - 2) + [F(-1)]
    return v[0] - v[-1]


def interval_gram(N, k):
    a = [[0] * N for _ in range(N)]
    for d in range(1, min(k, N) + 1):
        for start in range(N - d + 1):
            for i in range(start, start + d):
                for j in range(start, start + d):
                    a[i][j] += 1
    return a


def assert_psd(matrix):
    """Exact symmetric Schur-complement test, including zero pivots."""
    a = [[F(v) for v in row] for row in matrix]
    n = len(a)
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    for i in range(n):
        assert a[i][i] >= 0
        if not a[i][i]:
            assert all(a[i][j] == 0 for j in range(i, n))
            continue
        for j in range(i + 1, n):
            for h in range(j, n):
                a[j][h] -= a[i][j] * a[i][h] / a[i][i]
                a[h][j] = a[j][h]


def fib_lucas(m):
    f0, f1 = 0, 1
    l0, l1 = 2, 1
    for _ in range(m):
        f0, f1 = f1, f0 + f1
        l0, l1 = l1, l0 + l1
    return f0, l0


def fidelity():
    # arXiv:1712.05883v1, theorem labelled thm:effRes in the TeX:
    # m=n-2; r(1,n)=(m+1)/5 + 4 F_(m+1)/(5 L_(m+1)).
    for n in range(3, 31):
        f, l = fib_lucas(n - 1)
        expected = F(n - 1, 5) + F(4 * f, 5 * l)
        assert grounded_resistance(n, 2) == expected, (n, expected)
    print('FIDELITY: source Fibonacci/Lucas formula, k=2, n=3..30: 28 exact equalities')


def small_cases():
    total = defects = 0
    for k in range(1, 9):
        D = sum(d * d for d in range(1, k + 1))
        assert 6 * D == k * (k + 1) * (2 * k + 1)
        increments = []
        for n in range(2, 2 * k + 9):
            N = n - 1
            Q = interval_gram(N, k)
            x = solve(Q, [1] * N)
            R = grounded_resistance(n, k)
            assert sum(x, F(0)) == R, (n, k, 'electrical reduction')
            assert_psd([[Q[i][j] - int(i == j) for j in range(N)] for i in range(N)])
            assert_psd([[D * int(i == j) - Q[i][j] for j in range(N)] for i in range(N)])
            assert all(Q[i][j] == 0 for i in range(N) for j in range(N) if abs(i-j) > k-1)
            c = [D - sum(row) for row in Q]
            y = solve(Q, c)
            assert R == F(N, D) + (sum(c) + dot(c, y)) / D**2
            if N >= 2 * k:
                a = [sum(d * (d - i) for d in range(i + 1, k + 1))
                     if i < k else 0 for i in range(1, N + 1)]
                assert c == [a[i] + a[N - 1 - i] for i in range(N)]
                assert sum(v*v for v in c) == 2 * sum(v*v for v in a)
                defects += 1
            if k == 1:
                assert R == N
            increments.append(R)
            total += 1
        errors = [increments[i+1] - increments[i] - F(1, D)
                  for i in range(len(increments)-1)]
        print(f'SMALL_CASES: k={k}, n=2..{2*k+8}, cases={len(increments)}, '
              f'last_increment_error={errors[-1]}')
    print(f'IDENTITIES_AND_SPECTRAL_BOUNDS: {total} exact cases; stabilized defects: {defects}')


def boundary_moments():
    total = 0
    for k in range(2, 9):
        D = sum(d*d for d in range(1, k+1))
        b = k-1
        for ell in range(7):
            moments = []
            for N in ((ell+3)*b+2*k, (ell+3)*b+2*k+1):
                Q = interval_gram(N, k)
                B = [[D*int(i==j)-Q[i][j] for j in range(N)] for i in range(N)]
                a = [sum(d*(d-i) for d in range(i+1, k+1))
                     if i<k else 0 for i in range(1, N+1)]
                left, right = a[:], a[::-1]
                for _ in range(ell):
                    left, right = matvec(B, left), matvec(B, right)
                moments.append(dot(a, left))
                assert dot(a, right) == 0, (k, ell, N, 'opposite-boundary moment')
            assert moments[0] == moments[1], (k, ell, 'prefix stabilization')
            total += 1
    print(f'BOUNDARY_LOCALIZATION: {total} stabilized moments and 98 zero cross moments; '
          'k=2..8, powers=0..6; integer B=D*I-Q')


def main():
    fidelity()
    small_cases()
    boundary_moments()
    print('SCOPE: finite exact checks; no proof of the infinite limit or Lean theorem is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
