#!/usr/bin/env python3
"""Exact checks for arXiv:2611.06470v1 Section 7's three-point question.

Independently written from the paper's mathematics; no ancillary code used.
Finite-field ranks certify rational lower bounds. A dimension bound derived
from the S3 multiplicities supplies the matching upper bounds over Q/R.
The finite checks are evidence for the preregistered general argument, not a
proof of its universal quantifier.
"""
from collections import Counter
from itertools import permutations, product
from math import prod
import sympy as sp

PERMS = tuple(permutations(range(3)))
BASIS_ROWS = ((1, 0), (0, 1), (-1, -1))
P = tuple(tuple(BASIS_ROWS[p[i]] for i in range(2)) for p in PERMS)
EVEN, ODD = (0, 3, 4), (1, 2, 5)
PRIME = 65521


def excess(counts):
    return max(0, *(2 * counts[p] - sum(counts[q] for q in cls)
                    for cls in (EVEN, ODD) for p in cls))


def exact_rank_mod(rows):
    if not rows:
        return 0
    a = [[x % PRIME for x in row] for row in rows]
    r, n = 0, len(a[0])
    for c in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        u = pow(a[r][c], -1, PRIME)
        a[r] = [(x * u) % PRIME for x in a[r]]
        ar = a[r]
        for i in range(r + 1, len(a)):
            t = a[i][c]
            if t:
                a[i] = [(x - t * y) % PRIME for x, y in zip(a[i], ar)]
        r += 1
        if r == len(a):
            break
    return r


def full_matrix(labels, directions, k):
    return [[int(a == c) * d[i] for i in range(k) for c in range(2)] +
            [-P[p][a][c] * d[i] for i in range(k) for c in range(2)]
            for p, d in zip(labels, directions) for a in range(2)]


def axes(k):
    return [tuple(int(i == j) for i in range(k)) for j in range(k)]


def source_list_check(k, directions):
    """Axes eliminate X. Relabel X so the first permutation is identity.

    With delta the largest within-parity excess, on X=P Y the constraint
    count is at most 2k-delta; thus the reduced kernel has dimension >=delta.
    Rank >=2k-delta modulo PRIME therefore certifies equality over Q/R and
    containment in that global-permutation subspace. Delta=0 means kernel=0.
    """
    assert directions[:k] == axes(k)
    checked, fallback = 0, 0
    for tail in product(range(6), repeat=2*k-1):
        labels = (0,) + tail
        counts = Counter(labels)
        delta = excess(counts)
        rows = [[d[i] * (P[labels[i]][a][c] - P[p][a][c])
                 for i in range(k) for c in range(2)]
                for p, d in zip(labels[k:], directions[k:]) for a in range(2)]
        rank = exact_rank_mod(rows)
        target = 2*k - delta
        if rank != target:
            # A special rational source list need not be generic. Fall back to
            # the literal rational kernel-containment criterion from Section 6.
            fallback += 1
            m = sp.Matrix(rows)
            ker = m.nullspace()
            assert any(all(all(sum((P[labels[i]][a][c] - P[q][a][c]) * v[2*i+c]
                                  for c in range(2)) == 0
                               for i in range(k) for a in range(2)) for v in ker)
                       for q in range(6)), (k, labels, 'source list fails')
        checked += 1
    print(f'SOURCE k={k}: matching tuples={checked}; rational fallbacks={fallback}; separates triples')


def pair_reduce(c, cls):
    """Pair unlike elements within a parity; leave a dominant run or triad."""
    order = sorted(cls, key=lambda p: (-c[p], p))
    x, y, z = (c[p] for p in order)
    pairs = []
    residual = Counter()
    if x >= y + z:
        pairs += [(order[0], order[1])] * y
        pairs += [(order[0], order[2])] * z
        residual[order[0]] = x-y-z
    else:
        t = (x+y+z) % 2
        residual.update({p: t for p in cls})
        x, y, z = x-t, y-t, z-t
        for a, b, n in ((0, 1, (x+y-z)//2),
                        (0, 2, (x+z-y)//2),
                        (1, 2, (y+z-x)//2)):
            assert n >= 0
            pairs += [(order[a], order[b])] * n
    rec = Counter(residual)
    for a, b in pairs:
        assert a != b
        rec[a] += 1
        rec[b] += 1
    assert all(rec[p] == c[p] for p in cls)
    return pairs, +residual


def kind(c, cls):
    values = [c[p] for p in cls]
    support = sum(v > 0 for v in values)
    if support == 0:
        return 'empty', 0
    if support == 1:
        return 'run', sum(values)
    assert values == [1, 1, 1]
    return 'triad', 3


def classify(c):
    kinds = sorted((kind(c, EVEN), kind(c, ODD)))
    if all(t != 'triad' for t, _ in kinds):
        return 'A' if any(t == 'empty' for t, _ in kinds) else 'B'
    if all(t == 'triad' for t, _ in kinds):
        return 'E'
    run = next(n for t, n in kinds if t != 'triad')
    assert run % 2 == 1 and run >= 1
    return 'D' if run == 1 else 'C'


def compositions(n, length):
    if length == 1:
        yield (n,)
    else:
        for x in range(n+1):
            for tail in compositions(n-x, length-1):
                yield (x,) + tail


def reduction_check():
    seen, cases = 0, Counter()
    for k in range(0, 9):
        for counts in compositions(2*k, 6):
            c = Counter(dict(enumerate(counts)))
            pe, re = pair_reduce(c, EVEN)
            po, ro = pair_reduce(c, ODD)
            r = re + ro
            assert excess(c) == excess(r)
            assert 2*(len(pe)+len(po)) + sum(r.values()) == 2*k
            cases[classify(r)] += 1
            seen += 1
    print(f'REDUCTION all multiplicities for 0<=k<=8: {seen}; five cases {dict(sorted(cases.items()))}')
    assert set(cases) == set('ABCDE')


def matrix_checks():
    for cls in (EVEN, ODD):
        for p, q in permutations(cls, 2):
            assert abs((sp.Matrix(P[p])-sp.Matrix(P[q])).det()) == 3
    for p in EVEN:
        for q in ODD:
            assert (sp.Matrix(P[p])-sp.Matrix(P[q])).rank() == 1
    # A and B: a copies of I, b copies of one reflection, a>=b.
    for a in range(0, 25):
        for b in range(a+1):
            if (a+b) % 2:
                continue
            h = (a+b)//2
            if not h:
                continue
            es = axes(h)
            labels = [0]*a + [1]*b
            ds = es + [es[0]]*(a-h) + es[:b]
            assert exact_rank_mod(full_matrix(labels, ds, h)) == 4*h-a
    # C: odd run a>=3 opposite a triad; each of 3 constraints uses its own axis.
    for a in range(3, 26, 2):
        h = (a+3)//2
        es = axes(h)
        labels = [0]*a + list(ODD)
        ds = es + [es[0]]*(a-h) + es[:3]
        assert exact_rank_mod(full_matrix(labels, ds, h)) == 4*h-a
    # D: one I and three reflections; a 7x7 minor has determinant 3.
    d = sp.Matrix(full_matrix([0,1,2,5], [(1,0),(0,1),(1,1),(1,3)], 2))
    assert d[:7,:7].det() == 3
    assert d.rank() == 7
    assert all(v[:4, :] == v[4:, :] for v in d.nullspace())
    # E: all six permutations once; a 12x12 determinant is nonzero.
    e = sp.Matrix(full_matrix([0,3,4,1,2,5],
        axes(3)+[(1,2,4),(3,-1,2),(-2,5,1)], 3))
    assert e.det() == -19440
    # Structural alternatives: coordinate pairing for D and an edge-labelled
    # triangle for E. The ordinary proof uses ranks of permutation differences.
    d_axis = full_matrix([0,1,2,5], [(1,0),(1,0),(0,1),(0,1)], 2)
    assert exact_rank_mod(d_axis) == 7
    e_triangle = sp.Matrix(full_matrix([0,3,4,1,2,5],
        axes(3)+[(1,1,0),(0,1,1),(1,0,1)], 3))
    assert e_triangle.det() == 27
    print('RESIDUAL witnesses: A,B runs <=24; C odd runs <=25; D minor=3; E determinant=-19440')
    print('STRUCTURAL alternatives: D coordinate-pair rank=7; E triangle determinant=27')


def lower_bound_check():
    # Nonparallel branch of Theorem 1.2: two scaled orthogonal dipoles and
    # their sum. Coinciding middle terms cancel, leaving two disjoint triples.
    x = ((0,0),(2,1),(1,2))
    y = ((1,0),(0,1),(2,2))
    ds = ((0,1),(1,0),(1,-1))
    assert Counter(x) != Counter(y)
    for d in ds:
        assert sorted(sum(a*b for a,b in zip(v,d)) for v in x) == sorted(
            sum(a*b for a,b in zip(v,d)) for v in y)
    # The parallellism/zero-coefficient branches produce dipole products,
    # hence size <=2; pad both by the same point to obtain size exactly 3.
    print('LOWER source three-dipole cancellation: explicit distinct triples agree on three planar directions')


def main():
    assert sp.isprime(PRIME)
    lower_bound_check()
    matrix_checks()
    reduction_check()
    source_list_check(2, [(1,0),(0,1),(1,1),(1,3)])
    source_list_check(3, axes(3)+[(1,2,4),(3,-1,2),(-2,5,1)])
    source_list_check(4, axes(4)+[(5,2,-8,3),(-1,7,4,9),(8,-3,6,-5),(4,9,-2,7)])
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
