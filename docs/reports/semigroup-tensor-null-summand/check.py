#!/usr/bin/env python3
"""Exact checks for He--Tubbenhauer, arXiv:2508.04054v1, Conjecture 4.2.

Only Python integers are used. Binary matrices are tuples of bit-encoded rows.
The fidelity data are the four-vertex fusion matrix in Example 4.3.
Finite checks do not prove the universal conjecture.
"""

from collections import Counter
from itertools import product


def basis(vectors):
    pivots = {}
    for original in vectors:
        v = original
        while v:
            p = v.bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                break
    return tuple(pivots.values())


def rank(vectors):
    return len(basis(vectors))


def kernel(rows, dim):
    rows = list(rows)
    pivots = []
    h = 0
    for col in range(dim):
        pivot = next((i for i in range(h, len(rows)) if rows[i] >> col & 1), None)
        if pivot is None:
            continue
        rows[h], rows[pivot] = rows[pivot], rows[h]
        for i in range(len(rows)):
            if i != h and rows[i] >> col & 1:
                rows[i] ^= rows[h]
        pivots.append(col)
        h += 1
    return tuple((1 << f) | sum(((rows[i] >> f) & 1) << c
                               for i, c in enumerate(pivots))
                 for f in range(dim) if f not in pivots)


def transpose(a):
    return tuple(sum(((row >> j) & 1) << i for i, row in enumerate(a))
                 for j in range(len(a)))


def apply(a, v):
    return sum(((row & v).bit_count() & 1) << i for i, row in enumerate(a))


def compose(a, b):
    cols = tuple(apply(a, v) for v in transpose(b))
    return transpose(cols)


def kron(a, b):
    db = len(b)
    return tuple(sum((row_b if row_a >> j & 1 else 0) << (j * db)
                     for j in range(len(a))) for row_a in a for row_b in b)


def tensor_power(a, n):
    t = (1,)
    for _ in range(n):
        t = kron(t, a)
    return t


def member(v, span):
    return rank((*span, v)) == len(span)


def null_splitting(actions):
    """Compare kernel/image dimension with the dual pairing; construct the split."""
    dim = len(actions[0])
    k = kernel(tuple(row for a in actions for row in a), dim)
    image = basis(v for a in actions for v in transpose(a))
    intersection_dim = len(k) + len(image) - rank((*k, *image))
    mu = len(k) - intersection_dim
    dual = kernel(tuple(row for a in actions for row in transpose(a)), dim)
    pairing = tuple(sum(((f & v).bit_count() & 1) << i for i, v in enumerate(k))
                    for f in dual)
    assert rank(pairing) == mu

    # Select L inside K, disjoint from I. Extend I to U while remaining disjoint from L.
    l = []
    for v in k:
        if rank((*image, *l, v)) > len(image) + len(l):
            l.append(v)
    assert len(l) == mu
    u = list(image)
    for j in range(dim):
        v = 1 << j
        if rank((*l, *u, v)) > len(l) + len(u):
            u.append(v)
    assert len(l) + len(u) == dim and rank((*l, *u)) == dim
    assert all(apply(a, v) == 0 for a in actions for v in l)
    assert all(member(apply(a, v), u) for a in actions for v in u)
    assert len(image) <= sum(rank(a) for a in actions)
    assert dim - len(k) <= sum(rank(a) for a in actions)
    assert mu >= dim - 2 * sum(rank(a) for a in actions)
    return mu


def fidelity():
    # Column convention, ordered dimensions (1,6,1,8), third vertex the null module.
    a = ((0, 0, 0, 0), (1, 0, 0, 0), (0, 28, 6, 40), (0, 1, 0, 1))
    dims = (1, 6, 1, 8)
    assert all(sum(dims[i] * a[i][j] for i in range(4)) == 6 * dims[j]
               for j in range(4))
    v = (1, 0, 0, 0)
    counts = []
    for n in range(13):
        assert sum(d * m for d, m in zip(dims, v)) == 6 ** n
        counts.append(sum(v))
        if n >= 2:
            assert v == (0, 0, 6 ** n - 8, 1)
            assert sum(v) == 6 ** n - 7
        if n == 2:
            assert v[2] == 28
        v = tuple(sum(a[i][j] * v[j] for j in range(4)) for i in range(4))
    assert counts[:6] == [1, 1, 29, 209, 1289, 7769]
    print('FIDELITY Example 4.3: b(0..5)=', counts[:6])
    print('FIDELITY null multiplicity at n=2: 28; b(n)=6^n-7 for 2<=n<=12')


def binary_families():
    matrices = tuple(tuple(rows) for rows in product(range(4), repeat=2))
    singular = tuple(a for a in matrices if rank(a) < 2)
    assert len(singular) == 10
    zero = (0, 0)
    index = {a: i for i, a in enumerate(singular)}
    multiplication = tuple(tuple(index[compose(a, b)] for b in singular) for a in singular)
    tensors = {(i, n): tensor_power(a, n) for i, a in enumerate(singular)
               for n in (1, 2, 3)}
    for i, a in enumerate(singular):
        for n in (1, 2, 3):
            assert rank(tensors[i, n]) == rank(a) ** n
    tested = semigroups = nonunital = 0
    for mask in range(1, 1 << len(singular)):
        ids = tuple(i for i in range(len(singular)) if mask >> i & 1)
        actions = tuple(singular[i] for i in ids)
        closed = all(mask >> multiplication[i][j] & 1 for i in ids for j in ids)
        if closed:
            semigroups += 1
            has_identity = any(all(multiplication[e][i] == multiplication[i][e] == i
                                   for i in ids) for e in ids)
            nonunital += not has_identity
        # Independent maximal-summand decision: exhaust every binary 2x2 projection.
        allowed = tuple(e for e in matrices if compose(e, e) == e and
                        all(compose(a, e) == compose(e, a) == zero for a in actions))
        assert null_splitting(actions) == max(rank(e) for e in allowed)
        for n in (1, 2, 3):
            lifted = tuple(tensors[i, n] for i in ids)
            mu = null_splitting(lifted)
            r = sum(rank(singular[i]) ** n for i in ids)
            assert mu >= 2 ** n - 2 * r
            tested += 1
    print(f'EXHAUSTIVE GF(2): 1023 nonempty singular-matrix families; {tested} tensor cases n=1,2,3')
    print(f'Closed families: {semigroups}; closed families without identity: {nonunital}')
    print('All n=1 cases also match exhaustive enumeration of every null-image idempotent')
    # A uniform exact threshold for every subfamily of the ten singular matrices.
    n = 1
    while 2 ** n <= 2 * sum(rank(a) ** n for a in singular):
        n += 1
    mu = null_splitting(tuple(tensor_power(a, n) for a in singular))
    lower = 2 ** n - 2 * sum(rank(a) ** n for a in singular)
    assert lower > 0 and mu >= lower
    print(f'UNIFORM exact integer threshold n={n}: lower_bound={lower}, computed_mu={mu}')


def jordan_families():
    cases = 0
    for d in range(2, 7):
        for n in range(1, 6):
            vertices = tuple(product(range(d), repeat=n))
            successor = {x: tuple(i - 1 for i in x) for x in vertices if all(x)}
            incoming = set(successor.values())
            starts = set(vertices) - incoming
            blocks = Counter()
            visited = set()
            for start in starts:
                x = start
                length = 1
                visited.add(x)
                while x in successor:
                    x = successor[x]
                    assert x not in visited
                    visited.add(x)
                    length += 1
                blocks[length] += 1
            assert visited == set(vertices)
            assert sum(length * count for length, count in blocks.items()) == d ** n
            assert len(successor) == (d - 1) ** n
            assert sum(x in successor and successor[x] in successor for x in vertices) == (d - 2) ** n
            assert blocks[1] == d ** n - 2 * (d - 1) ** n + (d - 2) ** n
            assert sum(blocks.values()) == d ** n - (d - 1) ** n
            cases += 1
    print(f'JORDAN exact coordinate-chain decompositions: {cases} cases, 2<=d<=6, 1<=n<=5')
    print('m_Z=d^n-2(d-1)^n+(d-2)^n; b=d^n-(d-1)^n in every case')


if __name__ == '__main__':
    fidelity()
    binary_families()
    jordan_families()
    print('ALL_CHECKS_OK')
