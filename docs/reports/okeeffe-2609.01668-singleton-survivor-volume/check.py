#!/usr/bin/env python3
"""Exact checks for O'Keeffe arXiv:2609.01668v2, SM eq:sm:vsing.

Requires Python 3, sympy, and numpy. All evidence arithmetic is exact.
The Monte Carlo uses a finite rational grid; it is statistical corroboration,
not a proof of the continuous-volume identity.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import factorial, prod
import random
import numpy as np
import sympy as sp


def formula(n, gamma):
    return (1 - gamma) ** (n - 1) / (n * prod(n - j * gamma for j in range(1, n)))


def dot(u, v):
    return sum((a * b for a, b in zip(u, v)), Q(0))


def solve(a, b):
    """Rational elimination; None means a singular active-facet system."""
    n = len(b)
    m = [list(row) + [v] for row, v in zip(a, b)]
    for j in range(n):
        k = next((k for k in range(j, n) if m[k][j]), None)
        if k is None:
            return None
        m[j], m[k] = m[k], m[j]
        c = m[j][j]
        m[j] = [v / c for v in m[j]]
        for k in range(n):
            if k != j:
                c = m[k][j]
                m[k] = [x - c * y for x, y in zip(m[k], m[j])]
    return tuple(row[-1] for row in m)


def source_facets(n, a):
    """A*w <= b: ordered physical band and the first N-1 firing tests.

    The no-absorption test at a homogeneous state v is
    v[1] <= (1-q)*v[0]. Apply the source's L after each test.
    Equalities give the closure; they have zero Lebesgue measure.
    """
    m, q = 1 + a, a / (n * (1 + a))
    eye = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
    rows, rhs = [eye[0]], [m]
    for i in range(n - 1):
        rows.append([y - x for x, y in zip(eye[i], eye[i + 1])])
        rhs.append(Q(0))
    rows.append([-x for x in eye[-1]])
    rhs.append(Q(-1))
    power = eye
    for _ in range(n - 1):
        rows.append([y - (1 - q) * x for x, y in zip(power[0], power[1])])
        rhs.append(Q(0))
        power = [[power[i + 1][j] + q * power[0][j] for j in range(n)]
                 for i in range(n - 1)] + [[x / m for x in power[0]]]
    return rows, rhs


def transformed_vertices(n, a):
    """Inverse images of the N+1 weighted-simplex vertices."""
    m, s, q = 1 + a, a / n, a / (n * (1 + a))
    coeff = [1 - (n - 1 - j) * q for j in range(n)]
    assert all(c > 0 for c in coeff)
    out = []
    for k in range(-1, n):
        z = [Q(0)] * n
        if k >= 0:
            z[k] = s / coeff[k]
        h, ds = z[0], z[1:]
        y = [m - h - sum(ds[:j], Q(0)) for j in range(n)]
        w = tuple(y[j] - q * sum(y[:j], Q(0)) for j in range(n))
        out.append(w)
    return out


def simplex_probability(n, a, vertices):
    base = vertices[0]
    mat = sp.Matrix([[sp.Rational(vertices[j + 1][i] - base[i])
                      for j in range(n)] for i in range(n)])
    # simplex volume = |det|/N!, and the ordered w-density is N!/a^N.
    return Q(abs(mat.det())) / a ** n


def event_survival(x, gamma, events):
    """Direct rational voltage firing/reset dynamics, independent of L."""
    n = len(x)
    x = list(x)
    for _ in range(events):
        w = max(x)
        if x.count(w) != 1:
            return False
        j = x.index(w)
        y = [Q(0) if i == j else
             ((1 - gamma) * u + 1 - w) / (1 - gamma * w) + Q(1, n)
             for i, u in enumerate(x)]
        if any(v >= 1 for v in y):
            return False
        x = y
    return True


def exact_checks():
    # Fidelity: the source's explicit N=2 sync probability is (3-g)/(2*(2-g)).
    g = sp.symbols('g')
    assert sp.cancel(1 - (3 - g) / (2 * (2 - g)) -
                     (1 - g) / (2 * (2 - g))) == 0
    print('FIDELITY N=2: 1-P_sync=(1-gamma)/(2*(2-gamma)), symbolic residual=0')
    rng = random.Random(260901668)
    grids = 0
    for n in range(2, 11):
        for a in [Q(1, 4), Q(1, 2), Q(1), Q(4)]:
            rows, rhs = source_facets(n, a)
            expected = transformed_vertices(n, a)
            assert all(dot(row, v) <= b for v in expected
                       for row, b in zip(rows, rhs))
            p = simplex_probability(n, a, expected)
            assert p == formula(n, -a)
            # Independently enumerate the vertices of the original event polytope.
            if n <= 5 and a <= 1:
                actual = set()
                for inds in combinations(range(len(rows)), n):
                    v = solve([rows[i] for i in inds], [rhs[i] for i in inds])
                    if v is not None and all(dot(row, v) <= b
                                             for row, b in zip(rows, rhs)):
                        actual.add(v)
                assert actual == set(expected), (n, a, len(actual))
            # Unit-triangular recurrence, event inequalities and band identity.
            q, m, s = a / (n * (1 + a)), 1 + a, a / n
            for _ in range(30):
                w = sorted((Q(rng.randrange(1, 100000), 100000) * a + 1
                            for _ in range(n)), reverse=True)
                y, accum = [], Q(0)
                for v in w:
                    y.append(v + q * accum)
                    accum += y[-1]
                v = list(w)
                for k in range(n):
                    assert v[0] == y[k]
                    v = [v[i + 1] + q * v[0] for i in range(n - 1)] + [v[0] / m]
                assert v[0] == w[0] / m + q * sum(y[1:], Q(0))
                h = m - y[0]
                ds = [y[i - 1] - y[i] for i in range(1, n)]
                total = (1 - (n - 1) * q) * h + sum(
                    ((1 - (n - 1 - j) * q) * ds[j - 1]
                     for j in range(1, n)), Q(0))
                assert w[-1] - 1 == s - total
                x = [(v - 1) / a for v in w]
                predicate = all(y[i] > y[i + 1] for i in range(n - 1))
                assert predicate == event_survival(x, -a, n - 1)
                assert predicate == event_survival(x, -a, 3 * n)
                grids += 1
            # A guaranteed interior survivor, rather than only rare random survivors.
            coeff = [1 - (n - 1 - j) * q for j in range(n)]
            z = [s / (2 * n * c) for c in coeff]
            y = [m - z[0] - sum(z[1:j + 1], Q(0)) for j in range(n)]
            w = [y[j] - q * sum(y[:j], Q(0)) for j in range(n)]
            assert event_survival([(v - 1) / a for v in w], -a, 5 * n)
        # gamma=0: explicit staircase simplex, exact determinant and event dynamics.
        side = Q(1, n)
        base = [Q(n - 1 - i, n) for i in range(n)]
        vs = [base] + [[base[i] + (side if i < j else 0) for i in range(n)]
                       for j in range(1, n + 1)]
        determinant = sp.Matrix([[sp.Rational(vs[j + 1][i] - base[i])
                                 for j in range(n)] for i in range(n)]).det()
        assert abs(determinant) == Q(1, n ** n) == formula(n, Q(0))
        interior = [base[i] + Q(n - i, n * (n + 1)) for i in range(n)]
        assert event_survival(interior, Q(0), 5 * n)
        print(f'EXACT N={n}: 4 negative gamma values; gamma=0; '
              f'original-facet vertices independently enumerated={n <= 5}; OK')
    # The y-coordinate Jacobian is triangular with diagonal one, as is (h,d)->y.
    q = sp.symbols('q')
    for n in range(2, 11):
        ws = sp.symbols(f'w0:{n}')
        ys = []
        for w in ws:
            ys.append(w + q * sum(ys))
        jac = sp.Matrix(ys).jacobian(ws)
        assert jac.det() == 1
    print(f'EXACT direct-voltage/event-coordinate comparisons={grids}; no disagreements')
    print('EXACT unit-triangular Jacobian checks N=2..10: determinant=1')


def monte_carlo(n, trials=10000000):
    """gamma=-1/2; iid rational grid, exact integer homogeneous event maps.

    W=2*D+r encodes w=1+x/2 with x=r/D. The integer map is 3N*L:
    (TW)_i=3N*W_(i+1)+W_0, (TW)_last=2N*W_0.
    All products fit signed int64 through N-1 events (asserted below).
    The source's finite-firing criterion identifies survival statistically.
    """
    d = 2 ** 30
    assert 3 * d * (3 * n + 1) ** (n - 1) < 2 ** 63
    rng = np.random.Generator(np.random.PCG64(260901668 + n))
    hits = 0
    for start in range(0, trials, 100000):
        size = min(100000, trials - start)
        r = rng.integers(1, d, size=(size, n), dtype=np.int64)
        r.sort(axis=1)
        v = (2 * d + r[:, ::-1]).copy()
        # Initial equal voltages are a null-set analogue on the finite grid.
        distinct = np.all(v[:, :-1] != v[:, 1:], axis=1)
        v = v[distinct]
        for k in range(n - 1):
            v = v[(3 * n) * v[:, 1] < (3 * n - 1) * v[:, 0]]
            if not len(v):
                break
            if k + 1 < n - 1:
                v = np.column_stack(((3 * n) * v[:, 1:] + v[:, :1],
                                     2 * n * v[:, 0]))
        hits += len(v)
    p = formula(n, Q(-1, 2))
    z_squared = (Q(hits) - trials * p) ** 2 / (trials * p * (1 - p))
    # Fixed before sampling: flag a gross six-standard-deviation discrepancy.
    assert z_squared < 36, (n, hits, p, z_squared)
    print(f'MONTE_CARLO N={n} gamma=-1/2 seed={260901668+n} trials={trials} '
          f'grid_denominator={d} survivors={hits} observed={Q(hits,trials)} '
          f'formula={p} standardized_error_squared={z_squared} (<36)')


if __name__ == '__main__':
    exact_checks()
    monte_carlo(6)
    monte_carlo(7)
    # Include the nested 100M reading as well as the 10M pilot for the rarer N=7 event.
    monte_carlo(7, trials=100000000)
    print('ALL_CHECKS_OK')
