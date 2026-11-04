#!/usr/bin/env python3
"""Exact finite checks of the two OEIS A322844 matrix conjectures.

Sources: https://oeis.org/A322844 (revision 68),
https://oeis.org/A317614 (revision 194). No network is used by this check.
The checks supply finite evidence, not a proof for unbounded matrix order.
"""
from fractions import Fraction as Q
from itertools import combinations
import sys
import sympy as sp


def matrix(n):
    """Source's literal one-based entry formula, including the empty matrix."""
    return sp.Matrix(n, n, lambda i, j: (j + 1) + n * i if (i + 1) % 2
                     else n * (i + 1) - (j + 1) + 1)


def sequence(n):
    return Q(n * n * (3 * (1 + n * n) - 2 * (2 + n * n) * (n % 2)), 12)


def minor2(a, i, j):
    return int(a[i, i] * a[j, j] - a[i, j] * a[j, i])


def minor3(a, i, j, k):
    return int(a[i, i] * (a[j, j] * a[k, k] - a[j, k] * a[k, j])
               - a[i, j] * (a[j, i] * a[k, k] - a[j, k] * a[k, i])
               + a[i, k] * (a[j, i] * a[k, j] - a[j, j] * a[k, i]))


def main():
    print('Python ' + sys.version.split()[0] + '; SymPy ' + sp.__version__)
    assert matrix(2).tolist() == [[1, 2], [4, 3]]
    assert matrix(3).tolist() == [[1, 2, 3], [6, 5, 4], [7, 8, 9]]
    trace_terms = [1, 4, 15, 32, 65, 108, 175, 256, 369, 500]
    assert [int(sp.trace(matrix(n))) for n in range(1, 11)] == trace_terms
    scalar_terms = [0, 0, 5, 6, 68, 50, 333, 196, 1040, 540, 2525, 1210,
                    5220, 2366, 9653, 4200, 16448, 6936, 26325, 10830]
    assert [sequence(n) for n in range(20)] == scalar_terms
    print('FIDELITY: A317614 M(2), M(3), first 10 trace terms; A322844 first 20 terms: OK')
    totals = {'principal_2_minors': 0, 'principal_3_minors': 0,
              'higher_coefficients': 0, 'small_exhaustive_minors': 0}
    signed = []
    x = sp.Symbol('x')
    for n in range(41):
        a = matrix(n)
        epsilon = n % 2
        assert sp.trace(a) == Q(n ** 3 + n * epsilon, 2)
        if n > 2:
            assert a.det() == 0  # Fidelity: A317614 explicitly labels this proved.
        rank = a.rank()
        assert rank == min(n, 2)
        coeffs = a.charpoly(x).all_coeffs()  # Direct dense-matrix exact calculation.
        assert len(coeffs) == n + 1 and coeffs[0] == 1
        for k in range(3, n + 1):
            assert coeffs[k] == 0, (n, k, coeffs[k])
            totals['higher_coefficients'] += 1
        # k > n is the trace of an endomorphism of the zero exterior space.
        for k in range(3, 43):
            trace_ext = (-1) ** k * coeffs[k] if k <= n else 0
            assert trace_ext == 0, (n, k)
        pairs = list(combinations(range(n), 2))
        e2 = sum(minor2(a, i, j) for i, j in pairs)
        totals['principal_2_minors'] += len(pairs)
        assert e2 == -sequence(n)
        if n >= 2:
            assert coeffs[2] == e2 and abs(e2) == sequence(n)
            signed.append(e2)
        for i, j, k in combinations(range(n), 3):
            assert minor3(a, i, j, k) == 0, (n, i, j, k)
            totals['principal_3_minors'] += 1
        if n <= 10:
            for k in range(3, n + 1):
                for subset in combinations(range(n), k):
                    assert a.extract(subset, subset).det() == 0
                    totals['small_exhaustive_minors'] += 1
        # Independent rational checks of the proposed centered factorization and moments.
        w = [Q(i) - Q(n + 1, 2) for i in range(1, n + 1)]
        s = [(-1) ** (i + 1) for i in range(1, n + 1)]
        u = [n * i - Q(n - 1, 2) for i in range(1, n + 1)]
        assert all(Q(int(a[i, j])) == u[i] + s[i] * w[j]
                   for i in range(n) for j in range(n))
        assert sum(w) == 0
        assert sum(t * t for t in w) == Q(n * (n * n - 1), 12)
        assert sum(s) == epsilon
        assert sum(t * z for t, z in zip(w, s)) == -Q(n * (1 - epsilon), 2)
        vu = [[sum(u), sum(s)], [sum(t * z for t, z in zip(w, u)),
                                sum(t * z for t, z in zip(w, s))]]
        expected = [[Q(n * (n * n + 1), 2), epsilon],
                    [Q(n * n * (n * n - 1), 12), -Q(n * (1 - epsilon), 2)]]
        assert vu == expected
        assert vu[0][0] * vu[1][1] - vu[0][1] * vu[1][0] == e2
        print(f'n={n:2d} trace={int(sp.trace(a))} e2={e2} rank={rank} higher=0: OK')
    assert signed[:4] == [-5, -6, -68, -50]
    print('FIRST_SIGNED_E2: ' + ','.join(map(str, signed[:8])))
    print('COUNTS: ' + ', '.join(f'{k}={v}' for k, v in totals.items()))
    print('SCOPE: 0<=n<=40; all 3<=k<=42; direct principal 2- and 3-minors; all principal k>=3 minors for n<=10')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
