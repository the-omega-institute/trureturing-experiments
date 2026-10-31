#!/usr/bin/env python3
"""Exact source-fidelity and boundary checks for arXiv:2609.10403v1, Q9.1.

Theorem D supplies the independent three-vertex tree fidelity test.
The candidate check uses Theorems E/F's transfer/inverse-sparsity identities.
These are symbolic and rational computations, not a Lean settlement.
"""
from itertools import combinations, product
import sympy as s

Q = s.Rational
J = s.diag(1, -1, -1, -1, -1)
S = [0, 1, 2]
C = [0, 1, 2, 3]
D = [0, 1, 2, 4]
EDGES = [(i, j) for i, j in combinations(range(5), 2) if (i, j) != (3, 4)]


def eq(x, y):
    z = s.cancel(s.expand(x - y))
    assert not z.atoms(s.Float), z
    assert z == 0, z


def matrix_eq(A, B):
    assert A.shape == B.shape
    for x, y in zip(A, B):
        eq(x, y)


def gram(P):
    return P.T * J * P


def frame(t, e, a, b, ra, rb):
    return s.Matrix([
        [1, Q(5, 4), t, Q(5, 4), Q(5, 4)],
        [0, Q(3, 4), -Q(3, 4), 0, 0],
        [0, 0, e, a, b],
        [0, 0, 0, ra, 0],
        [0, 0, 0, 0, rb],
    ])


def principal(A, indices):
    return A.extract(indices, indices)


def fidelity():
    # Theorem D: the product-distance completion on the path 1--2--3.
    # All three points are future unit vectors in R^(1,2).
    U = s.Matrix([[Q(5, 4), 1, Q(13, 12)],
                  [Q(3, 4), 0, 0], [0, 0, Q(5, 12)]])
    B = U.T * s.diag(1, -1, -1) * U
    assert list(B.diagonal()) == [1, 1, 1]
    assert all(U[0, j] > 0 for j in range(3))
    eq(B[0, 1], Q(5, 4))
    eq(B[1, 2], Q(13, 12))
    eq(B[0, 2], B[0, 1] * B[1, 2])
    # Proposition 2.7(1): determinant is the product of (1-a_e^2).
    eq(B.det(), (1 - B[0, 1] ** 2) * (1 - B[1, 2] ** 2))
    assert B.det() != 0
    eq(B.inv()[0, 2], 0)
    print('FIDELITY: Theorem D path product 65/48; determinant 25/256; inverse nonedge zero')


def graph_check():
    # The only maximal cliques are C and D. An explicit perfect elimination
    # ordering is also checked. Chordality is independently decided by
    # testing every cyclic ordering on every subset of size at least four.
    V = set(range(5))
    adj = {frozenset(e) for e in EDGES}
    def clique(X):
        return all(frozenset(e) in adj for e in combinations(X, 2))
    cliques = [set(X) for k in range(1, 6) for X in combinations(range(5), k) if clique(X)]
    maximal = [X for X in cliques if not any(X < Y for Y in cliques)]
    assert {frozenset(X) for X in maximal} == {frozenset(C), frozenset(D)}
    remaining = V.copy()
    for v in [3, 4, 0, 1, 2]:
        assert clique([w for w in remaining if frozenset([v, w]) in adj])
        remaining.remove(v)
    from itertools import permutations
    for k in range(4, 6):
        for c in permutations(range(5), k):
            cycle = {frozenset([c[i], c[(i + 1) % k]]) for i in range(k)}
            assert {e for e in adj if e <= set(c)} != cycle
    assert len(EDGES) == 9
    assert all(frozenset([0, j]) in adj for j in range(1, 5))
    print('GRAPH: connected, 9 edges, exactly two maximal cliques, no induced cycle of length >=4')


def symbolic_family():
    t, e, a, b, ra, rb, x = s.symbols('t e a b ra rb x', real=True)
    P = frame(t, e, a, b, ra, rb)
    A = gram(P)
    # Relations specify t=sqrt(25/16+e^2), ra=sqrt(9/16-a^2),
    # rb=sqrt(9/16-b^2), with the positive roots used in actual instances.
    def normalized(z):
        return s.expand(z).subs({t ** 2: Q(25, 16) + e ** 2,
                                ra ** 2: Q(9, 16) - a ** 2,
                                rb ** 2: Q(9, 16) - b ** 2}).expand()
    for v in A.diagonal():
        eq(normalized(v), 1)
    eq(P.det(), Q(3, 4) * e * ra * rb)
    sep = principal(A, S)
    eq(sep.det(), Q(9, 16) * e ** 2)
    eq(principal(A, C).det(), -Q(9, 16) * e ** 2 * ra ** 2)
    eq(principal(A, D).det(), -Q(9, 16) * e ** 2 * rb ** 2)
    T = P.extract(S, S)
    matrix_eq(sep, T.T * s.diag(1, -1, -1) * T)
    wa, wb = P.extract(S, [3]), P.extract(S, [4])
    matrix_eq(A.extract(S, [3]), T.T * s.diag(1, -1, -1) * wa)
    matrix_eq(A.extract(S, [4]), T.T * s.diag(1, -1, -1) * wb)
    # Theorem E, with one clique-tree edge.
    transfer = (A.extract([3], S) * sep.inv() * A.extract(S, [4]))[0]
    eq(transfer, Q(25, 16) - a * b)
    Ax = A.copy()
    Ax[3, 4] = Ax[4, 3] = x
    residual = Ax.extract([3, 4], [3, 4]) - Ax.extract([3, 4], S) * sep.inv() * Ax.extract(S, [3, 4])
    matrix_eq(residual, s.Matrix([[-ra ** 2, x - Q(25, 16) + a * b],
                                 [x - Q(25, 16) + a * b, -rb ** 2]]))
    eq(Ax.det(), Q(9, 16) * e ** 2 * (ra ** 2 * rb ** 2 - (x - Q(25, 16) + a * b) ** 2))
    # Theorem F's inverse-sparsity characterization, at the transfer value.
    eq(A.inv()[3, 4], 0)
    eq(Q(9, 16) - a ** 2, (Q(3, 4) - a) * (Q(3, 4) + a))
    print('SYMBOLIC: unit relations, frame/clique/separator determinants, transfer 25/16-a*b, Schur certificate, inverse sparsity')


def rational_sequences():
    m = s.symbols('m', integer=True, positive=True)
    t = 5 * (m ** 2 + 1) / (4 * (m ** 2 - 1))
    e = 5 * m / (2 * (m ** 2 - 1))
    eq(t ** 2 - e ** 2, Q(25, 16))
    eq(t - Q(5, 4), 5 / (2 * (m ** 2 - 1)))
    assert s.limit(t, m, s.oo) == Q(5, 4)
    assert s.limit(e, m, s.oo) == 0
    # Both denominators are positive for m>=2; these exact identities also
    # certify e>0 and t-e=5(m-1)/(4(m+1))>0 on that domain.
    eq(t - e, 5 * (m - 1) / (4 * (m + 1)))
    eq(Q(5, 4) * t - Q(9, 20) * e - Q(4, 5) * t, Q(9, 20) * (t - e))
    A0 = s.Matrix([
        [1, Q(5, 4), Q(5, 4), Q(5, 4), Q(5, 4)],
        [Q(5, 4), 1, Q(17, 8), Q(25, 16), Q(25, 16)],
        [Q(5, 4), Q(17, 8), 1, Q(25, 16), Q(25, 16)],
        [Q(5, 4), Q(25, 16), Q(25, 16), 1, 0],
        [Q(5, 4), Q(25, 16), Q(25, 16), 0, 1],
    ])  # The two zeros at positions (4,5),(5,4) are unspecified, not data.
    limits = []
    for a, b, ra, rb in [(Q(0), Q(0), Q(3, 4), Q(3, 4)),
                         (Q(9, 20), Q(9, 20), Q(3, 5), Q(3, 5))]:
        A = gram(frame(t, e, a, b, ra, rb))
        L = A.applyfunc(lambda z: s.limit(z, m, s.oo))
        limits.append(L)
        for i in range(5):
            eq(L[i, i], 1)
        for i, j in EDGES:
            eq(L[i, j], A0[i, j])
            assert L[i, j] > 1
        P0 = frame(Q(5, 4), Q(0), a, b, ra, rb)
        matrix_eq(L, gram(P0))
        assert all(P0[0, j] > 0 for j in range(5))
        assert len({tuple(P0[:, j]) for j in range(5)}) == 5
        matrix_eq(P0[:, 2], Q(5, 2) * P0[:, 0] - P0[:, 1])
        assert principal(L, C).rank() == principal(L, D).rank() == 3
        eq(principal(L, [0, 1, 3]).det(), Q(81, 256))
        eq(principal(L, [0, 1, 4]).det(), Q(81, 256))
        eq(A[3, 4], Q(25, 16) - a * b)
    eq(limits[0][3, 4], Q(25, 16))
    eq(limits[1][3, 4], Q(34, 25))
    eq(limits[0][3, 4] - limits[1][3, 4], Q(81, 400))
    print('LIMITS: exact rational functions; same prescribed A0; clique ranks 3; five distinct limiting points; all prescribed entries >1')
    print('COUNTEREXAMPLE: canonical limits 25/16 and 34/25; difference 81/400')


def finite_crosscheck():
    # Rational points on a circle of radius 3/4. These checks are independent
    # evaluations of the symbolic certificate, not evidence for an unbounded
    # universal obtained by enumeration.
    us = [Q(-3, 4), Q(-1, 2), Q(-1, 3), Q(0), Q(1, 3), Q(1, 2), Q(3, 4)]
    tested = 0
    for m, u, v in product(range(2, 7), us, us):
        t, e = Q(5 * (m * m + 1), 4 * (m * m - 1)), Q(5 * m, 2 * (m * m - 1))
        a, b = Q(3, 4) * 2 * u / (1 + u * u), Q(3, 4) * 2 * v / (1 + v * v)
        ra, rb = Q(3, 4) * (1 - u * u) / (1 + u * u), Q(3, 4) * (1 - v * v) / (1 + v * v)
        P = frame(t, e, a, b, ra, rb)
        A = gram(P)
        assert all(P[0, j] > 0 for j in range(5))
        assert list(A.diagonal()) == [1] * 5
        assert ra > 0 and rb > 0 and e > 0
        assert all(A[i, j] > 1 for i, j in EDGES)
        assert P.det() != 0 and A.det() > 0
        assert principal(A, C).det() < 0 and principal(A, D).det() < 0
        transfer = (A.extract([3], S) * principal(A, S).inv() * A.extract(S, [4]))[0]
        eq(transfer, Q(25, 16) - a * b)
        eq(A[3, 4], transfer)
        eq(A.inv()[3, 4], 0)
        assert not A.atoms(s.Float)
        tested += 1
    print(f'EXACT_CROSSCHECK: {tested} rational instances, m=2..6 and 7x7 circle parameters')


if __name__ == '__main__':
    fidelity()
    graph_check()
    symbolic_family()
    rational_sequences()
    finite_crosscheck()
    print('ALL_CHECKS_OK')
