#!/usr/bin/env python3
"""Exact checks for Snellman, arXiv:2609.23605v1, Conjecture 52.

All numerical calculations use integers or fractions.Fraction. No numerical LP,
floating point hull routine, or third-party source code is used.
"""
from fractions import Fraction as F
from itertools import combinations


def dot(a, x):
    return sum((p * q for p, q in zip(a, x)), F(0))


def rref(rows, ncols):
    a = [[F(v) for v in row] for row in rows]
    pivots = []
    for j in range(ncols):
        k = next((k for k in range(len(pivots), len(a)) if a[k][j]), None)
        if k is None:
            continue
        i = len(pivots)
        a[i], a[k] = a[k], a[i]
        c = a[i][j]
        a[i] = [v / c for v in a[i]]
        for k in range(len(a)):
            if k != i and a[k][j]:
                c = a[k][j]
                a[k] = [v - c * w for v, w in zip(a[k], a[i])]
        pivots.append(j)
    return a, pivots


def rank(rows, ncols):
    return len(rref(rows, ncols)[1])


def solve(rows, rhs):
    n = len(rows[0])
    a, pivots = rref([list(r) + [b] for r, b in zip(rows, rhs)], n)
    if len(pivots) != n:
        return None
    assert all(not any(row[:n]) or not row[n] for row in a[n:])
    return tuple(a[i][n] for i in range(n))


def basis(n, i):
    return tuple(F(j == i) for j in range(n))


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def affine_dimension(points, n):
    points = list(points)
    return rank([sub(x, points[0]) for x in points[1:]], n) if points else -1


def staircase(n, d, i, shifted=False):
    return tuple(F((d if shifted else 0) <= j < (d if shifted else 0) + i)
                 for j in range(n + d))


def halfspaces(n, d):
    """Return Ax >= b in the source x coordinates; x_1=1 is separate."""
    N = n + d
    result = []
    for j in range(N - 1):
        row = sub(basis(N, j), basis(N, j + 1))
        result.append((f'chain_{j + 1}', row, F(-1 if j + 1 == d else 0)))
    result.append(('floor', basis(N, N - 1), F(0)))
    result.append(('entry_cap', basis(N, d), F(1)))
    result.append(('exit_cap', tuple(-v for v in basis(N, n)), F(-1)))
    return result


def feasible(x, hs):
    return x[0] == 1 and all(dot(a, x) >= b for _, a, b in hs)


def enumerate_vertices(n, d, hs):
    """Every full-rank choice of N-1 active inequalities, over Q."""
    N = n + d
    result = set()
    trials = 0
    for active in combinations(hs, N - 1):
        trials += 1
        x = solve([basis(N, 0)] + [a for _, a, _ in active],
                  [F(1)] + [b for _, _, b in active])
        if x is not None and feasible(x, hs):
            result.add(x)
    assert result
    return result, trials


def difference(x):
    return tuple(x[j] - x[j + 1] if j + 1 < len(x) else x[j]
                 for j in range(len(x)))


def transform(x, d):
    return add(difference(x), basis(len(x), d - 1))


def inverse(z, d):
    w = sub(z, basis(len(z), d - 1))
    return tuple(sum(w[j:], F(0)) for j in range(len(z)))


def block_masses(z, n, d):
    return sum(z[:d], F(0)), sum(z[d:n], F(0)), sum(z[n:], F(0))


def split_mass(z, n, d):
    a, b, c = block_masses(z, n, d)
    assert all(v >= 0 for v in z) and a + b + c == 2
    assert a <= 1 and c <= 1
    if b:
        theta = (1 - a) / b
        assert 0 <= theta <= 1
        u = z[:d] + tuple(theta * q for q in z[d:n]) + (F(0),) * d
        v = (F(0),) * d + tuple((1 - theta) * q for q in z[d:n]) + z[n:]
    else:
        assert a == c == 1
        u = z[:d] + (F(0),) * n
        v = (F(0),) * n + z[n:]
    assert all(q >= 0 for q in u + v)
    assert sum(u) == sum(v) == 1 and add(u, v) == z
    assert all(u[j] == 0 for j in range(n, n + d))
    assert all(v[j] == 0 for j in range(d))
    return u, v


def check_affine_identities(n, d, hs):
    N = n + d
    # An affine expression on x_1=1 is determined by these N test points.
    anchor = basis(N, 0)
    tests = [anchor] + [add(anchor, basis(N, j)) for j in range(1, N)]
    for x in tests:
        z = transform(x, d)
        assert inverse(z, d) == x
        a, b, c = block_masses(z, n, d)
        assert a + b + c == x[0] + 1 == 2
        assert a == x[0] - x[d] + 1 and c == x[n]
        slacks = [dot(row, x) - rhs for _, row, rhs in hs]
        assert slacks == list(z) + [1 - a, 1 - c]
    # Thus H maps into {z>=0, sum z=2}, proving boundedness of H.
    # The following identities use the input vertices independently of H.
    for i in range(1, n + 1):
        vi = staircase(n, d, i)
        ui = staircase(n, d, i, shifted=True)
        assert difference(vi) == basis(N, i - 1)
        assert difference(ui) == sub(basis(N, d + i - 1), basis(N, d - 1))


def check_case(n, d):
    N, m = n + d, n - d
    hs = halfspaces(n, d)
    sums = {add(staircase(n, d, i), staircase(n, d, j, shifted=True))
            for i in range(1, n + 1) for j in range(1, n + 1)}
    assert all(feasible(x, hs) for x in sums)
    vertices, trials = enumerate_vertices(n, d, hs)
    # Bounded H and these two inclusions give H = conv(all pair sums)
    # = conv(unshifted vertices) + conv(shifted vertices).
    assert vertices <= sums
    assert affine_dimension(vertices, N) == N - 1
    assert len(vertices) == n * n - m * (m - 1)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            x = add(staircase(n, d, i), staircase(n, d, j, shifted=True))
            t = i - d  # signed integer, as in the source theorem
            predicted = not (1 <= t <= m and 1 <= j <= m and t != j)
            assert (x in vertices) == predicted
    facets = {}
    face_dimensions = {}
    for name, a, b in hs:
        face = frozenset(x for x in vertices if dot(a, x) == b)
        dim = affine_dimension(face, N)
        face_dimensions[name] = dim
        if dim == N - 2:
            assert face not in facets
            facets[face] = name
    assert len(facets) == (N + 1 if m == 1 else N + 2)
    if m == 1:
        assert face_dimensions[f'chain_{d + 1}'] == N - 3
    else:
        assert all(dim == N - 2 for dim in face_dimensions.values())
    check_affine_identities(n, d, hs)
    # Exact mass-split checks include interior points and the b=0 boundary.
    samples = list(vertices)
    samples.append(tuple(sum(x[j] for x in vertices) / len(vertices) for j in range(N)))
    samples.extend(tuple((x[j] + y[j]) / 2 for j in range(N))
                   for x, y in combinations(sorted(vertices), 2))
    boundary = add(basis(N, 0), basis(N, n))
    assert block_masses(boundary, n, d) == (1, 0, 1)
    split_mass(boundary, n, d)
    for x in samples:
        assert feasible(x, hs)
        split_mass(transform(x, d), n, d)
    if (n, d) == (5, 2):
        assert (affine_dimension(vertices, N), len(vertices), len(facets)) == (6, 19, 9)
        assert all(x[0] == 1 and x[0] != 5 for x in sums)
        print('FIDELITY Example 53: dimension=6 vertices=19 facets=9')
        print('SOURCE_TYPO Example 53 x_1=5 conflicts with Definition 36; x_1=1 checked')
        print('FIDELITY Theorem 51: extremal matrix and vertex formula checked')
    print(f'n={n} d={d} N={N} active_bases={trials} '
          f'vertices={len(vertices)} facets={len(facets)} '
          f'mass_split_samples={len(samples) + 1} H_EQUALS_MINKOWSKI_OK')
    return trials, len(samples) + 1


def main():
    cases = bases = samples = 0
    for n in range(2, 8):
        for d in range(1, n):
            count, splits = check_case(n, d)
            cases += 1
            bases += count
            samples += splits
    print(f'EXACT_SUMMARY cases={cases} active_bases={bases} mass_splits={samples}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
