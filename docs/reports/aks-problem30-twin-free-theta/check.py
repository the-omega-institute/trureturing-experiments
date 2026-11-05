#!/usr/bin/env python3
"""Exact pre-probe checks for Alfuraidan--Khan--Sauras-Altuzarra Problem 30.

All matrices come directly from edge lists. Characteristic polynomials and
nullities use exact integer/rational arithmetic. No numerical eigensolver is used.
"""
from collections import deque
from itertools import combinations
import sympy as sp

X = sp.Symbol('x')


def graph(vertices, edges):
    vertices = tuple(vertices)
    assert len(vertices) == len(set(vertices))
    pos = {v: i for i, v in enumerate(vertices)}
    neighbors = [set() for _ in vertices]
    for u, v in edges:
        i, j = pos[u], pos[v]
        assert i != j and j not in neighbors[i]
        neighbors[i].add(j)
        neighbors[j].add(i)
    adjacency = sp.zeros(len(vertices))
    for i, ns in enumerate(neighbors):
        for j in ns:
            adjacency[i, j] = 1
    degree = sp.diag(*(len(ns) for ns in neighbors))
    distances = []
    for start in range(len(vertices)):
        row = [None] * len(vertices)
        row[start] = 0
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in neighbors[i]:
                if row[j] is None:
                    row[j] = row[i] + 1
                    queue.append(j)
        assert all(d is not None for d in row), 'graph must be connected'
        distances.append(row)
    return vertices, pos, neighbors, distances, degree - adjacency, degree + adjacency


def theta(lengths):
    vertices = ['u', 'v']
    edges = []
    for i, length in enumerate(lengths):
        assert isinstance(length, int) and length >= 2
        interior = [(i, j) for j in range(1, length)]
        vertices.extend(interior)
        path = ['u', *interior, 'v']
        edges.extend(zip(path, path[1:]))
    return graph(vertices, edges)


def twin_free(data):
    neighbors = data[2]
    return all(neighbors[i] - {j} != neighbors[j] - {i}
               for i, j in combinations(range(len(neighbors)), 2))


def resolves(data, landmarks):
    distances = data[3]
    profiles = [tuple(distances[v][s] for s in landmarks)
                for v in range(len(distances))]
    return len(set(profiles)) == len(profiles)


def metric_dimension(data, upper=3):
    for size in range(upper + 1):
        for landmarks in combinations(range(len(data[0])), size):
            if resolves(data, landmarks):
                return size
    raise AssertionError('no resolving set within the declared bound')


def spectral(data):
    values, polynomials = [], []
    for matrix in data[4:6]:
        assert matrix == matrix.T
        poly = matrix.charpoly(X).as_poly()
        assert poly.domain == sp.ZZ
        # In characteristic zero each squarefree factor has distinct roots.
        # Symmetry makes all roots real and geometric = algebraic multiplicity.
        _, factors = sp.sqf_list(poly)
        values.append(max(exponent for _, exponent in factors))
        polynomials.append(poly)
    return tuple(values), polynomials


def fidelity_checks():
    # The source's Table 2 (label TwinFreeExamples), rows 1 and 4.
    rows = [
        (1, [(1, 2), (1, 6), (1, 7), (2, 3), (2, 6),
             (3, 4), (4, 5), (5, 6)], (2, 2, 1)),
        (4, [(1, 2), (1, 5), (1, 6), (2, 3), (2, 6),
             (2, 7), (3, 4), (3, 6), (4, 5)], (2, 1, 2)),
    ]
    for row, edges, expected in rows:
        data = graph(range(1, 8), edges)
        mm, _ = spectral(data)
        actual = (metric_dimension(data), *mm)
        assert twin_free(data) and actual == expected
        print(f'FIDELITY twin-free table row={row}: (dim,mmL,mmQ)={actual}')
    # Theorem 25, Theta_{2,2,r}, including both congruence directions.
    for r in (3, 5, 7, 9):
        data = theta((2, 2, r))
        mm, _ = spectral(data)
        expected = (2, 1) if r % 4 == 3 else (1, 2)
        assert metric_dimension(data) == 2 and mm == expected
        assert not twin_free(data)
        print(f'FIDELITY Theorem 25 r={r}: dim=2, mmL={mm[0]}, mmQ={mm[1]}, twins=True')


def symmetry_restrictions(data, r):
    vertices, pos = data[:2]
    representatives = ['u', 'v', (0, 1), (0, 2)] + [(2, j) for j in range(1, r)]
    plus = sp.zeros(len(vertices), len(representatives))
    for col, v in enumerate(representatives):
        plus[pos[v], col] = 1
        if v in ((0, 1), (0, 2)):
            plus[pos[(1, v[1])], col] = 1
    minus = sp.zeros(len(vertices), 2)
    for j in (1, 2):
        minus[pos[(0, j)], j - 1] = 1
        minus[pos[(1, j)], j - 1] = -1
    assert plus.row_join(minus).rank() == len(vertices)
    restrictions = []
    for matrix, sign in zip(data[4:6], (-1, 1)):
        product = matrix * plus
        w = product.extract([pos[v] for v in representatives], range(plus.cols))
        anti = sp.Matrix([[2, sign], [sign, 2]])
        assert matrix * plus == plus * w
        assert matrix * minus == minus * anti
        assert sp.gcd(w.charpoly(X).as_poly(), w.charpoly(X).as_poly().diff()).degree() == 0
        restrictions.append(w)
    return restrictions


def family_checks():
    for k in range(1, 9):
        r = 6 * k + 4
        h = r // 2
        data = theta((3, 3, r))
        vertices, pos, neighbors, distances = data[:4]
        assert len(vertices) == r + 5 == 6 * k + 9
        assert sum(map(len, neighbors)) // 2 == r + 6
        assert min(map(len, neighbors)) == 2
        assert len(neighbors[pos['u']]) == len(neighbors[pos['v']]) == 3
        assert twin_free(data)
        landmarks = [pos[(0, 1)], pos[(2, 1)]]
        assert resolves(data, landmarks)
        assert not resolves(data, ())
        assert all(not resolves(data, (i,)) for i in range(len(vertices)))
        short = ['u', 'v', (0, 1), (0, 2), (1, 1), (1, 2)]
        expected_short = [(1, 1), (2, 4), (0, 2), (1, 3), (2, 2), (3, 3)]
        assert [tuple(distances[pos[v]][s] for s in landmarks) for v in short] == expected_short
        for j in range(1, r):
            profile = tuple(distances[pos[(2, j)]][s] for s in landmarks)
            assert profile == (min(j + 1, r - j + 2), min(j - 1, r - j + 4))
            if j <= h:
                expected = (j + 1, j - 1)
            elif j == h + 1:
                expected = (h + 1, h)
            elif j == h + 2:
                expected = (h, h + 1)
            else:
                expected = (r - j + 2, r - j + 4)
            assert profile == expected
        mm, (chi_l, chi_q) = spectral(data)
        repeated_l = sp.gcd(chi_l, chi_l.diff()).monic()
        repeated_q = sp.gcd(chi_q, chi_q.diff()).monic()
        assert sp.expand(repeated_l.as_expr() - (X - 1) * (X - 3)) == 0
        assert repeated_q.as_expr() == 1
        assert mm == (2, 1)
        w_l, w_q = symmetry_restrictions(data, r)
        for eigenvalue in (1, 3):
            assert len((data[4] - eigenvalue * sp.eye(len(vertices))).nullspace()) == 2
            assert len((data[5] - eigenvalue * sp.eye(len(vertices))).nullspace()) == 1
            assert len((w_l - eigenvalue * sp.eye(w_l.rows)).nullspace()) == 1
            assert len((w_q - eigenvalue * sp.eye(w_q.rows)).nullspace()) == 0
        print(f'FAMILY k={k}, r={r}, n={len(vertices)}: twin_free=True, dim=2, mmL=2, mmQ=1; gcdL={repeated_l.as_expr()}, gcdQ=1; W simple, W_L mult(1,3)=(1,1), W_Q mult(1,3)=(0,0)')


def symbolic_checks():
    t = sp.Symbol('t')
    p = [sp.Poly(1, t), sp.Poly(t, t)]
    for _ in range(2, 81):
        p.append(sp.Poly(t * p[-1].as_expr() - p[-2].as_expr(), t))
    for n in range(1, len(p)):
        a, b = p[n].as_expr(), p[n - 1].as_expr()
        assert sp.expand(a*a + b*b - t*a*b - 1) == 0
        assert sp.expand(a.subs(t, -t) - (-1)**n * a) == 0
        assert p[n].eval(2) == n + 1
        assert p[n].eval(-2) == (-1)**n * (n + 1)
    for epsilon in (1, -1):
        denominator = t*t - 1
        forced_p = -denominator / 2
        forced_r = -((t + epsilon)*denominator - 2*t) / 2
        obstruction = sp.factor(forced_p**2 + forced_r**2 - t*forced_p*forced_r - 1)
        if epsilon == 1:
            expected = (t - 2)*(t - 1)*(t + 1)*(t*t + 2*t - 1)/4
            quadratic, value = t*t + 2*t - 1, 1
        else:
            expected = -(t - 1)*(t + 1)*(t + 2)*(t*t - 2*t - 1)/4
            quadratic, value = t*t - 2*t - 1, -1
        assert sp.expand(obstruction - expected) == 0
        assert sp.Poly(quadratic, t).is_irreducible
        assert sp.rem(forced_r - value, quadratic, t) == 0
        print(f'SYMBOLIC epsilon={epsilon}: obstruction={obstruction}; forced R={value} modulo {quadratic}')
    assert sp.Matrix([[1, -1], [1, 0]])**6 == sp.eye(2)
    assert sp.Matrix([[-1, -1], [1, 0]])**3 == sp.eye(2)
    print('SYMBOLIC continuant identities n=1..80; parity; values at +/-2; singular periods 6 and 3: OK')


def main():
    print(f'Python exact arithmetic; SymPy {sp.__version__}; no floating-point evidence')
    fidelity_checks()
    family_checks()
    symbolic_checks()
    print('Finite checks support the proposed infinite-family argument; they do not prove its universal quantifier.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
