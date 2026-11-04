#!/usr/bin/env python3
"""Exact exhaustive check of Cazenave arXiv:2611.03879v1, Conjecture 5.

Requires Python 3, SymPy, and nauty 2.8.9 geng. No floating point is used.
Every connected graph up to isomorphism of orders 2 through 9 is generated.
Integer symmetric congruence computes rank(A) and inertia((2m-rank(A))I-L).
At most one negative eigenvalue of the latter is exactly the target inequality.
"""
import argparse
import hashlib
from itertools import product
from math import gcd
from pathlib import Path
import subprocess
import sympy as sp

COUNTS = {2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117, 9: 261080}


def inertia(matrix):
    """Return positive, negative, zero inertia by exact integer congruence.

    A nonzero diagonal is a 1x1 pivot. If all diagonals vanish, a nonzero
    off-diagonal gives a 2x2 pivot with one sign of each kind. The Schur
    complement is multiplied by the absolute pivot, a strictly positive
    number. Dividing the entire remainder by its positive gcd controls size.
    """
    a = [row[:] for row in matrix]
    pos = neg = 0
    while a:
        size = len(a)
        k = next((i for i in range(size) if a[i][i]), None)
        if k is not None:
            pivot = a[k][k]
            sign = 1 if pivot > 0 else -1
            pos += sign > 0
            neg += sign < 0
            rest = [i for i in range(size) if i != k]
            b = [[0] * len(rest) for _ in rest]
            for ii, i in enumerate(rest):
                for jj in range(ii, len(rest)):
                    j = rest[jj]
                    v = abs(pivot) * a[i][j] - sign * a[i][k] * a[k][j]
                    b[ii][jj] = b[jj][ii] = v
        else:
            pair = next(((i, j) for i in range(size) for j in range(i + 1, size)
                         if a[i][j]), None)
            if pair is None:
                return pos, neg, size
            k, l = pair
            pivot = a[k][l]
            sign = 1 if pivot > 0 else -1
            pos += 1
            neg += 1
            rest = [i for i in range(size) if i not in pair]
            b = [[0] * len(rest) for _ in rest]
            for ii, i in enumerate(rest):
                for jj in range(ii, len(rest)):
                    j = rest[jj]
                    v = abs(pivot) * a[i][j] - sign * (
                        a[i][k] * a[l][j] + a[i][l] * a[k][j])
                    b[ii][jj] = b[jj][ii] = v
        common = 0
        for row in b:
            for v in row:
                common = gcd(common, v)
        if common > 1:
            b = [[v // common for v in row] for row in b]
        a = b
    return pos, neg, 0


def sympy_inertia(a):
    """Independent characteristic-polynomial rational root isolation."""
    p = n = z = 0
    for (lo, hi), mult in sp.Matrix(a).charpoly().as_poly().intervals(
            eps=sp.Rational(1, 1024)):
        if lo == hi == 0:
            z += mult
        elif lo >= 0:
            p += mult
        elif hi <= 0:
            n += mult
        else:
            raise AssertionError('root isolation must separate a nonzero root from zero')
    assert p + n + z == len(a)
    return p, n, z


def decode_graph6(line):
    raw = line.strip()
    n = raw[0] - 63
    assert 2 <= n <= 9
    values = [x - 63 for x in raw[1:]]
    assert all(0 <= x < 64 for x in values)
    bits = [(x >> shift) & 1 for x in values for shift in range(5, -1, -1)]
    assert len(bits) >= n * (n - 1) // 2
    a = [[0] * n for _ in range(n)]
    t = 0
    for j in range(1, n):
        for i in range(j):
            a[i][j] = a[j][i] = bits[t]
            t += 1
    assert not any(bits[t:])
    reached = {0}
    todo = [0]
    while todo:
        i = todo.pop()
        for j in range(n):
            if a[i][j] and j not in reached:
                reached.add(j)
                todo.append(j)
    assert len(reached) == n
    return a


def laplacian(a):
    return [[sum(a[i]) if i == j else -a[i][j]
             for j in range(len(a))] for i in range(len(a))]


def check_graph(a):
    n = len(a)
    p, q, nullity = inertia(a)
    rank = p + q
    assert rank + nullity == n
    m = sum(map(sum, a)) // 2
    bound = 2 * m - rank
    l = laplacian(a)
    shifted = [[(bound if i == j else 0) - l[i][j]
                for j in range(n)] for i in range(n)]
    signs = inertia(shifted)
    assert sum(signs) == n
    assert signs[1] <= 1, (a, rank, m, signs)
    equal = signs[1] + signs[2] >= 2
    return rank, m, equal, shifted


def fidelity():
    graphs = {
        'K2': [[0, 1], [1, 0]],
        'K3': [[0, 1, 1], [1, 0, 1], [1, 1, 0]],
        'P3': [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
        'P4': [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]],
        'K1,3': [[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]],
    }
    expected = {'K2': (0, 2, 1), 'K3': (3, 3, 3), 'P3': (1, 2, 2),
                'P4': (2, 4, 3), 'K1,3': (1, 2, 3)}
    for name, a in graphs.items():
        l = sp.Matrix(laplacian(a))
        roots = sorted([v for v, mult in l.eigenvals().items() for _ in range(mult)])
        mu2 = roots[-2]
        rank, m, equal, _ = check_graph(a)
        assert (mu2, rank, m) == expected[name]
        assert rank == sp.Matrix(a).rank()
        slack = 2 * m - rank - mu2
        assert equal == (slack == 0)
        print(f'FIDELITY {name}: mu2={mu2}; rank={rank}; m={m}; slack={slack}; '
              f'lap_charpoly={sp.factor(l.charpoly().as_expr())}', flush=True)
    print('SOURCE_EQUALITIES_OK K2 K3 P4', flush=True)
    for values in product((-1, 0, 1), repeat=6):
        x, y, z, u, v, w = values
        a = [[x, u, v], [u, y, w], [v, w, z]]
        assert inertia(a) == sympy_inertia(a), a
    print('INERTIA_CROSSCHECK_OK 729 symmetric integer 3x3 matrices', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geng', default='geng', help='nauty 2.8.9 geng executable')
    args = parser.parse_args()
    print('ARITHMETIC exact integers; SymPy exact algebra and rational root isolation', flush=True)
    print('CHECK_SHA256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), flush=True)
    fidelity()
    total = equality_count = 0
    stream_hash = hashlib.sha256()
    for n, expected in COUNTS.items():
        proc = subprocess.Popen([args.geng, '-cq', str(n)], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        count = eq = 0
        seen = set()
        for raw in proc.stdout:
            assert raw not in seen
            seen.add(raw)
            stream_hash.update(raw)
            a = decode_graph6(raw)
            assert len(a) == n
            rank, m, equal, shifted = check_graph(a)
            if n <= 5:
                assert rank == sp.Matrix(a).rank()
                assert inertia(shifted) == sympy_inertia(shifted)
            if equal:
                eq += 1
                print(f'EQUALITY n={n}; graph6={raw.strip().decode()}; rank={rank}; m={m}',
                      flush=True)
            count += 1
        error = proc.stderr.read().decode()
        assert proc.wait() == 0, error
        assert count == expected, (n, count, expected)
        total += count
        equality_count += eq
        print(f'ORDER n={n}; connected_unlabelled={count}; verified={count}; '
              f'failures=0; zero_slack={eq}', flush=True)
    assert total == 273192
    assert equality_count == 3
    print('GRAPH6_STREAM_SHA256', stream_hash.hexdigest(), flush=True)
    print(f'EXHAUSTIVE_OK connected simple graphs 2<=n<=9; total={total}; failures=0', flush=True)
    print('MINIMUM_SLACK 0 (exact; attained by K2, K3 and P4)', flush=True)
    print('BOUNDARY no conclusion for n>9 from this computation; no Lean proof', flush=True)
    print('ALL_CHECKS_OK', flush=True)


if __name__ == '__main__':
    main()
