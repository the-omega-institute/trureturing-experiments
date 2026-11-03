#!/usr/bin/env python3
"""Exact independent checks of Lin's definitions and a sharpness obstruction.

Requires Python 3 and SymPy. All evidence uses integers or exact rationals;
randomness chooses rational inputs only. No floating point arithmetic is used.
"""
from itertools import combinations
from random import Random
import sympy as S


def inertia(M):
    """Exact symmetric congruence elimination, including zero diagonal pivots."""
    M = S.Matrix(M)
    assert M == M.T
    positive = negative = 0
    while M.rows:
        diag = next((i for i in range(M.rows) if M[i, i] != 0), None)
        if diag is not None:
            order = [diag] + [i for i in range(M.rows) if i != diag]
            M = M.extract(order, order)
            pivot = M[0, 0]
            positive += int(bool(pivot > 0))
            negative += int(bool(pivot < 0))
            M = M[1:, 1:] - M[1:, :1] * M[:1, 1:] / pivot
        else:
            pair = next(((i, j) for i in range(M.rows)
                         for j in range(i + 1, M.rows) if M[i, j] != 0), None)
            if pair is None:
                return positive, negative, M.rows
            i, j = pair
            order = [i, j] + [k for k in range(M.rows) if k not in pair]
            M = M.extract(order, order)
            P = M[:2, :2]
            assert P.det() < 0
            positive += 1
            negative += 1
            M = M[2:, 2:] - M[2:, :2] * P.inv() * M[:2, 2:]
    return positive, negative, 0


def connected(vertices, edges):
    if not vertices:
        return False
    reached = {next(iter(vertices))}
    while True:
        new = reached | {v for u, v in edges if u in reached and v in vertices}
        if new == reached:
            return reached == set(vertices)
        reached = new


def powerset(xs):
    xs = list(xs)
    for size in range(len(xs) + 1):
        for cs in combinations(xs, size):
            yield frozenset(cs)


def weak_covers(A, y):
    """Exhaust Definition 5.9, including omitted zeros and zero-only members.

    Contract homogeneous/zero-zero edges first. Distinct zero domains have no
    edge. Conditions 1 and 4 then force each member's edges from its vertex set:
    all zero-strong edges, and no strong-strong edges. Enumerate these subgraphs,
    assigning each strong domain exactly once and testing all five conditions.
    Count domain combinations as sets of inverse-image vertex sets.
    """
    n = A.rows
    edges = {(i, j) for i in range(n) for j in range(n)
             if i != j and A[i, j] != 0}
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for i, j in edges:
        if y[i] * A[i, j] * y[j] < 0 or y[i] == y[j] == 0:
            parent[root(i)] = root(j)
    domains = tuple(frozenset(i for i in range(n) if root(i) == k)
                    for k in sorted({root(i) for i in range(n)}))
    phi = {v: d for d, vs in enumerate(domains) for v in vs}
    vertices = frozenset(range(len(domains)))
    zeros = frozenset(d for d in vertices if all(y[v] == 0 for v in domains[d]))
    strong = vertices - zeros
    H = {(phi[i], phi[j]) for i, j in edges if phi[i] != phi[j]}
    assert not any(i in zeros and j in zeros for i, j in H)
    allowed = {(i, j) for i, j in H if i in zeros or j in zeros}
    members = [vs for vs in powerset(vertices) if connected(vs, allowed)]
    by_strong = {s: [vs for vs in members if s in vs] for s in strong}
    accepted = set()

    def neighbor(vs, z):
        return frozenset(t for t in vs if (z, t) in allowed)

    def nonzero_sign(v, target):
        terms = [A[v, u] * y[u] for d in target for u in domains[d]
                 if (v, u) in edges]
        return ((any(t > 0 for t in terms) and all(t >= 0 for t in terms)) or
                (any(t < 0 for t in terms) and all(t <= 0 for t in terms)))

    def valid(C):
        # Clause 1; induced allowed edges were forced by clauses 1 and 4.
        if not all(connected(vs, allowed) for vs in C):
            return False
        # Clause 2.
        if not all(sum(s in vs for vs in C) == 1 for s in strong):
            return False
        # Clause 3, including exact cardinalities and sign unanimity.
        retained = set(strong)
        for z in zeros:
            at_z = [vs for vs in C if z in vs]
            N = frozenset(t for t in vertices if (z, t) in H)
            if not at_z:
                continue
            retained.add(z)
            if len(at_z) == 1:
                if neighbor(at_z[0], z) != N:
                    return False
            elif len(at_z) == 2:
                X, Y = [neighbor(vs, z) for vs in at_z]
                if X & Y or X | Y != N:
                    return False
                if not any(nonzero_sign(v, X) or nonzero_sign(v, Y)
                           for v in domains[z]):
                    return False
            else:
                return False
        # Clause 4.
        if not all((s, z) in allowed for vs in C for s in vs & strong
                   for z in vs & zeros if (s, z) in H):
            return False
        # Clause 5 uses the original contracted graph, not the union of members.
        return connected(retained, H)

    def finish(C):
        # Enumerate every possible zero-only connected member: a singleton.
        for singles in powerset(zeros):
            candidate = frozenset(C) | {frozenset([z]) for z in singles}
            if valid(candidate):
                combination = frozenset(frozenset().union(*(domains[d] for d in vs))
                                        for vs in candidate)
                accepted.add(combination)

    def assign(C, assigned):
        if assigned == strong:
            finish(C)
            return
        s = min(strong - assigned)
        for vs in by_strong[s]:
            if vs & assigned:
                continue
            if any(sum(z in w for w in C) + int(z in vs) > 2 for z in zeros):
                continue
            assign(C + [vs], assigned | (vs & strong))

    assign([], frozenset())
    assert accepted
    return domains, accepted


def family(r):
    m = r + 1
    E = S.Matrix(r, m, lambda i, j: r if i == j else -1)
    A = S.BlockMatrix([[S.zeros(r), E], [E.T, S.zeros(m)]]).as_explicit()
    y = S.Matrix([0] * r + [1] * m)
    return E, A, y


def main():
    # Fidelity: Lin Example 5.3 and Example 5.11 (label ex:wdone).
    A = S.zeros(7)
    A[0, 0] = 6
    for i in range(1, 7):
        A[0, i] = A[i, 0] = 1
    for i, j in [(1, 2), (3, 4), (5, 6)]:
        A[i, j] = A[j, i] = 1
    y = S.Matrix([0, 1, -1, 1, -1, 1, -1])
    assert A * y == -y
    assert A.eigenvals() == {-1: 3, 0: 1, 1: 2, 7: 1}
    domains, covers = weak_covers(A, y)
    assert len(domains) == 4 and len(covers) == 1
    assert max(map(len, covers)) == 1
    assert inertia(A + S.eye(7)) == (4, 0, 3)
    print('FIDELITY Lin Examples 5.3/5.11: spectrum={-1:3,0:1,1:2,7:1}; Ay=-y; W=1; n_-(A+I)+1=1')

    E, A, y = family(2)
    assert A * y == S.zeros(5, 1)
    assert A.eigenvals() == {-3: 1, -S.sqrt(3): 1, 0: 1, S.sqrt(3): 1, 3: 1}
    assert inertia(A) == (2, 2, 1)
    domains, covers = weak_covers(A, y)
    assert all(len(d) == 1 for d in domains)
    assert max(map(len, covers)) == 2
    assert {len(c) for c in covers} == {1, 2}
    print('WITNESS r=2: A=', A.tolist())
    print('WITNESS r=2: spectrum={-3,-sqrt(3),0,sqrt(3),3}; inertia=(2,2,1); '
          f'weak_combinations={len(covers)}; W=2; n_-+1=3')

    t = S.Symbol('t')
    for r in range(2, 8):
        E, A, y = family(r)
        assert E * E.T == (r + 1)**2 * S.eye(r) - (r + 1) * S.ones(r)
        assert A * y == S.zeros(2 * r + 1, 1)
        expected = t * (t*t - (r+1)) * (t*t - (r+1)**2)**(r-1)
        assert S.expand(A.charpoly(t).as_expr() - expected) == 0
        assert inertia(A) == (r, r, 1)
        if r <= 3:
            _, covers = weak_covers(A, y)
            assert max(map(len, covers)) == 2
            print(f'FAMILY r={r}: exact weak_combinations={len(covers)}; W=2; gap={r-1}')
    print('FAMILY r=2..7: exact Gram identity, characteristic polynomial and inertia all match')

    rng = Random(261100789)
    count = 0
    for r in range(2, 8):
        for sample in range(12):
            v = S.Matrix([S.Rational(rng.randrange(1, 10), rng.randrange(1, 8))
                          for _ in range(r + 1)])
            F = S.Matrix(r, r + 1, lambda i, j:
                         0 if i == j else -S.Rational(rng.randrange(1, 10), rng.randrange(1, 8)))
            for i in range(r):
                F[i, i] = -sum(F[i, j] for j in range(r + 1) if j != i)
            B = F * S.diag(*(1 / x for x in v))
            D = S.diag(*(S.Rational(rng.randrange(-9, 10), rng.randrange(1, 8))
                         for _ in range(r)))
            M = S.BlockMatrix([[D, B], [B.T, S.zeros(r+1)]]).as_explicit()
            assert B * v == S.zeros(r, 1)
            assert all(F[i, i] - sum(abs(F[i, j]) for j in range(r) if j != i)
                       == -F[i, r] > 0 for i in range(r))
            assert B.rank() == r
            assert inertia(M) == (r, r, 1)
            original = S.BlockMatrix([[S.zeros(r), B], [B.T, S.zeros(r+1)]]).as_explicit()
            assert inertia(original) == (r, r, 1)
            _, template, _ = family(r)
            assert all(S.sign(original[i, j]) == S.sign(template[i, j])
                       for i in range(2*r+1) for j in range(2*r+1))
            # Independently compute an explicit congruence to r hyperbolic planes.
            P = B[:, :r]
            Q = P.inv().row_join(-P.inv() * B[:, r:])
            Q = Q.col_join(S.zeros(1, r).row_join(S.eye(1)))
            assert Q.det() != 0 and B * Q == S.eye(r).row_join(S.zeros(r, 1))
            T = S.diag(S.eye(r), Q)
            N = T.T * M * T
            U = S.eye(2*r+1)
            U[r:2*r, :r] = -D / 2
            target = S.zeros(2*r+1)
            target[:r, r:2*r] = target[r:2*r, :r] = S.eye(r)
            assert U.T * N * U == target
            if r == 2 and sample < 4:
                _, covers = weak_covers(original, S.zeros(r, 1).col_join(v))
                assert max(map(len, covers)) == 2
            count += 1
    print(f'RANDOM seed=261100789: {count} original-pattern and {count} diagonal-relaxed exact rational realizations, r=2..7; '
          'Bv=0, strict dominance, full row rank, inertia=(r,r,1), explicit hyperbolic congruence')
    print('RANDOM first four r=2 realizations: independently enumerated W=2')
    print('LIMIT: finite exact checks support the proposed general argument; no Lean settlement proof is included')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
