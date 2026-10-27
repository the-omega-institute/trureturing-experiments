#!/usr/bin/env python3
"""Exact source fidelity and finite counterexample checks; Python standard library."""
from collections import Counter, deque
from functools import reduce
from itertools import product
from operator import xor

# Batal, arXiv:2609.05934v1, Theorem 4.3 neighborhood table.
SOURCE27 = [
    [7,9,11,14,18,25], [8,11,12,13,17,23], [5,6,10,19,20,24],
    [7,8,9,15,19,26], [3,6,10,16,22,27], [3,5,12,14,23,27],
    [1,4,13,19,21,26], [2,4,15,18,21,24], [1,4,16,20,21,22],
    [3,5,12,13,25,26], [1,2,16,17,22,26], [2,6,10,17,24,25],
    [2,7,10,15,23,27], [1,6,15,17,20,25], [4,8,13,14,18,27],
    [5,9,11,18,20,21], [2,11,12,14,22,23], [1,8,15,16,20,24],
    [3,4,7,21,22,27], [3,9,14,16,18,24], [7,8,9,16,19,25],
    [5,9,11,17,19,23], [2,6,13,17,22,24], [3,8,12,18,20,23],
    [1,10,12,14,21,26], [4,7,10,11,25,27], [5,6,13,15,19,26],
]
# Complete candidate open neighborhoods; vertices are one-based.
CANDIDATE16 = [
    [3,6,9,10,11,12], [4,5,8,9,10,15], [1,7,8,10,12,13],
    [2,6,7,9,13,14], [2,6,7,8,12,16], [1,4,5,11,14,16],
    [3,4,5,11,13,15], [2,3,5,11,14,16], [1,2,4,13,15,16],
    [1,2,3,14,15,16], [1,6,7,8,14,15], [1,3,5,13,14,15],
    [3,4,7,9,12,16], [4,6,8,10,11,12], [2,7,9,10,11,12],
    [5,6,8,9,10,13],
]
TREE = [(1,3),(1,6),(1,9),(1,10),(1,11),(1,12),(3,7),(3,8),
        (3,13),(6,4),(6,5),(6,14),(6,16),(9,2),(9,15)]
# (surviving variables, indices of D_i to sum), D_i = N(i) XOR N(1).
CERTIFICATES = [
    ((1,10,15,16), (3,4,6,8,11,13)),
    ((2,10), (2,6,8,9,10,12,14)),
    ((3,10,15,16), (2,6,7,12,14)),
    ((4,10), (2,3,6,7,8,15)),
    ((5,10,15,16), (2,4,9,11)),
    ((6,10), (4,5,8,11,12,15)),
    ((7,10,15,16), (4,5,7,9,10,11,12)),
    ((8,16), (10,11,15)),
    ((9,15), (3,8,10,11,12)),
    ((11,15), (3,4,11,15)),
    ((12,16), (4,5,9,10,14)),
    ((13,15), (4,6,8,11,14)),
    ((14,16), (2,3,12,13)),
]
EXPECTED_SCALARS = [0,639,21845,22314,43221,43690,64896,65535]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mask(indices):
    return sum(1 << (i-1) for i in indices)


def xor_sum(values):
    return reduce(xor, values, 0)


def graph_checks(table, degree, total):
    n = len(table)
    require(all(len(vs) == len(set(vs)) == degree for vs in table), 'degree/duplicates')
    require(all(1 <= j <= n and j != i for i, vs in enumerate(table, 1) for j in vs),
            'range/loops')
    require(all((j in table[i-1]) == (i in table[j-1])
                for i in range(1, n+1) for j in range(1, n+1)), 'asymmetry')
    require([sum(vs) for vs in table] == [total]*n, 'ordinary magic sums')
    # Independent reconstruction from unordered edges.
    edges = [(i,j) for i, vs in enumerate(table,1) for j in vs if i < j]
    rebuilt = [set() for _ in table]
    for i,j in edges:
        rebuilt[i-1].add(j)
        rebuilt[j-1].add(i)
    require(rebuilt == list(map(set, table)), 'edge reconstruction')
    seen = {1}
    queue = deque([1])
    while queue:
        for j in rebuilt[queue.popleft()-1]:
            if j not in seen:
                seen.add(j)
                queue.append(j)
    require(len(seen) == n, 'disconnected')
    return edges


def rank_mod(matrix, p):
    a = [[v % p for v in row] for row in matrix]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [(v*inv) % p for v in a[r]]
        for i in range(len(a)):
            if i != r:
                k = a[i][c]
                a[i] = [(x-k*y) % p for x,y in zip(a[i],a[r])]
        r += 1
        if r == len(a):
            break
    return r


def determinant_bareiss(matrix):
    a = [row[:] for row in matrix]
    previous, sign = 1, 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k,len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        v = a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                numerator = a[i][j]*v-a[i][k]*a[k][j]
                require(numerator % previous == 0, 'nonexact Bareiss division')
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = v
    return sign*a[-1][-1]


def labels(a,b,c):
    d = a ^ b ^ c
    return [d,a,d,a,d,a,d,c,b,a,b,c,b,c,b,c]


def generated_subgroup(values):
    # Additive closure in (Z/2Z)^4, represented by four-bit integers.
    subgroup = {0}
    for v in values:
        subgroup |= {x ^ v for x in subgroup}
    require(all(x ^ y in subgroup for x in subgroup for y in subgroup), 'closure')
    return subgroup


def main():
    source_edges = graph_checks(SOURCE27, 6, 84)
    a = [[int(j in vs) for j in range(1,28)] for vs in SOURCE27]
    reduced = [[a[i][j]-a[26][j] for j in range(26)] for i in range(26)]
    minor = [row[:24] for row in reduced[2:]]
    require(determinant_bareiss(minor) == -11, 'source minor determinant')
    require(rank_mod(reduced,3) == 24, 'source reduced ternary rank')
    print('FIDELITY Theorem 4.3: vertices=27 edges=%d degree=6 connected=True' % len(source_edges))
    print('FIDELITY identity-label neighborhood sums: [84] * 27')
    print('FIDELITY reduced F_3 rank=24 nullity=2; 24x24 minor determinant=-11')
    edges = graph_checks(CANDIDATE16,6,51)
    require(len(edges) == 48, 'edge count')
    require(len(TREE) == 15 and len(set(TREE)) == 15, 'tree edge count')
    require(all(tuple(sorted(e)) in edges for e in TREE), 'tree edge membership')
    tree_seen = {1}
    while True:
        enlarged = tree_seen | {v for u,v in TREE if u in tree_seen} | {u for u,v in TREE if v in tree_seen}
        if enlarged == tree_seen:
            break
        tree_seen = enlarged
    require(tree_seen == set(range(1,17)), 'spanning tree')
    print('COUNTEREXAMPLE vertices=16 edges=48 degree=6 connected=True; spanning tree checked')
    print('COUNTEREXAMPLE identity-label neighborhood sums: [51] * 16')
    rows = list(map(mask,CANDIDATE16))
    differences = [row ^ rows[0] for row in rows]
    for target, indices in CERTIFICATES:
        require(xor_sum(differences[i-1] for i in indices) == mask(target), 'certificate cancellation')
    require(rank_mod([[(row >> j)&1 for j in range(16)] for row in differences],2) == 13,
            'binary difference rank')
    print('CERTIFICATES all 13 XOR cancellations checked; difference rank=13')
    scalars = [x for x in range(1 << 16)
               if len({(x & row).bit_count() % 2 for row in rows}) == 1]
    formula_scalars = sorted(mask(i+1 for i,v in enumerate(labels(a,b,c)) if v)
                             for a,b,c in product(range(2),repeat=3))
    require(scalars == EXPECTED_SCALARS == formula_scalars, 'exhaustive scalar classification')
    print('SCALAR all 65536 binary maps tested; magic masks=%s' % scalars)
    histogram = Counter()
    for a,b,c in product(range(16),repeat=3):
        xs = labels(a,b,c)
        require({xor_sum(xs[j-1] for j in vs) for vs in CANDIDATE16} == {a ^ b},
                'vector formula magic sums')
        generated = generated_subgroup(xs)
        subset_sums = {xor_sum(v for v,bit in zip((a,b,c),bits) if bit)
                       for bits in product(range(2),repeat=3)}
        require(generated == subset_sums and len(generated) <= 8, 'three-generator ceiling')
        histogram[len(generated)] += 1
    require(dict(histogram) == {1:1,2:105,4:1470,8:2520}, 'subgroup histogram')
    require(len({tuple(labels(a,b,c)) for a,b,c in product(range(16),repeat=3)}) == 4096,
            'distinct vector formula maps')
    print('VECTOR all 4096 parameter triples checked; common sum=a XOR b')
    print('VECTOR generated subgroup size histogram=%s; maximum=8 < target order 16' % dict(sorted(histogram.items())))
    print('SCOPE exact integer/binary computations; no Lean proof; no minimal-order or infinite-family claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
