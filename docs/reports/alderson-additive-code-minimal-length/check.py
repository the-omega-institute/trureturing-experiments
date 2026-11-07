#!/usr/bin/env python3
"""Exact finite checks for Alderson, arXiv:2607.22297v1, Remark 7.8.

Original checker; Python standard library only. F_2^6 vectors are six-bit integers,
and F_2^3 alphabet symbols are three-bit integers. Addition is XOR. No floating
point operations, external certificates, or Lean proof are used.
"""

from collections import Counter
from itertools import combinations, permutations


def span(generators):
    result = {0}
    for g in generators:
        result |= {x ^ g for x in result}
    return frozenset(result)


def matchings(points):
    if not points:
        yield ()
        return
    a = points[0]
    for b in points[1:]:
        for rest in matchings(tuple(x for x in points if x not in (a, b))):
            yield ((a, b),) + rest


def rho(matching, epsilon, x):
    a, b = zip(*matching)
    pa = [(x >> i) & 1 for i in a]
    pb = [(x >> i) & 1 for i in b]
    s = [pa[i] ^ pb[i] for i in range(3)]
    return (s[0] ^ s[1]) | ((s[0] ^ s[2]) << 1) | (
        (pb[0] ^ pb[1] ^ pb[2] ^ (epsilon & s[0])) << 2
    )


def null_plane(matching, epsilon):
    ds = [(1 << a) ^ (1 << b) for a, b in matching]
    t = sum(1 << a for a, _ in matching)
    return span((ds[0] ^ ds[1], ds[0] ^ ds[2], t ^ (epsilon * ds[0])))


def code(coordinates):
    words = [tuple(rho(m, e, x) for m, e in coordinates) for x in range(64)]
    assert len(set(words)) == 64
    assert all({w[j] for w in words} == set(range(8)) for j in range(len(coordinates)))
    assert all(
        words[x ^ y] == tuple(a ^ b for a, b in zip(words[x], words[y]))
        for x in range(64) for y in range(64)
    )
    return words


def distances(words):
    return Counter(sum(a != b for a, b in zip(x, y)) for x, y in combinations(words, 2))


def extension(words, tile_plane):
    tile = {0, 63} | {1 << i for i in range(6)}
    labels = {}
    for label, w in enumerate(sorted(tile_plane)):
        for t in tile:
            assert t ^ w not in labels
            labels[t ^ w] = label
    assert set(labels) == set(range(64))
    assert {t ^ u for t in tile for u in tile}.isdisjoint(F)
    # The alphabet labels on the plane respect its F_2 addition.
    assert all(labels[w ^ v] == labels[w] ^ labels[v] for w in tile_plane for v in tile_plane)
    assert any(labels[x ^ y] != labels[x] ^ labels[y] for x in range(64) for y in range(64))
    return [word + (labels[x],) for x, word in enumerate(words)]


V = range(64)
F = frozenset(x for x in V if x.bit_count() in (3, 4))
assert len(F) == 35
all_spaces = {span(gs) for gs in combinations(range(1, 64), 3)}
three_spaces = {w for w in all_spaces if len(w) == 8}
assert len(three_spaces) == 1395
assert all(w & F for w in three_spaces)
planes = {w for w in three_spaces if w - {0} <= F}
assert len(planes) == 30
assert all(Counter(x.bit_count() for x in w - {0}) == {3: 4, 4: 3} for w in planes)
assert Counter(sum(x in w for w in planes) for x in F) == {6: 35}

source_coordinates = [(m, e) for m in matchings(tuple(range(6))) for e in (0, 1)]
assert len(source_coordinates) == 30
assert {null_plane(m, e) for m, e in source_coordinates} == planes
for m, e in source_coordinates:
    assert {x for x in V if rho(m, e, x) == 0} == set(null_plane(m, e))
    assert {rho(m, e, x) for x in V} == set(range(8))
source = code(source_coordinates)
source_histogram = distances(source)
assert source_histogram == {24: 1120, 30: 896}
source_tile_plane = span((int('110100', 2), int('011010', 2), int('001101', 2)))
assert source_tile_plane in planes
source_extension_histogram = distances(extension(source, source_tile_plane))
assert source_extension_histogram == {25: 1120, 30: 224, 31: 672}
print('FIDELITY: 30 source null planes; 35 forbidden points, each covered 6 times')
print('FIDELITY: source code cardinality=64; every coordinate has 8 symbols')
print('FIDELITY: source pair distances=' + str(dict(sorted(source_histogram.items()))))
print('FIDELITY: source extension pair distances=' + str(dict(sorted(source_extension_histogram.items()))))

five = [((5, i), ((i + 1) % 5, (i - 1) % 5), ((i + 2) % 5, (i - 2) % 5)) for i in range(5)]
two_coordinates = [(m, e) for m in five for e in (0, 1)]
two_cover = [null_plane(m, e) for m, e in two_coordinates]
assert len(set(two_cover)) == 10
assert all(sum(x in w for w in two_cover) == (2 if x in F else 0) for x in range(1, 64))
short = code(two_coordinates)
short_histogram = distances(short)
assert short_histogram == {8: 1120, 10: 896}
short_extension_histogram = distances(extension(short, source_tile_plane))
assert short_extension_histogram == {9: 1120, 10: 224, 11: 672}
print('REFUTATION CHECK: faithful F_2-linear (10,2,8)_(8/2) code, cardinality=64')
print('REFUTATION CHECK: pair distances=' + str(dict(sorted(short_histogram.items()))))
print('REFUTATION CHECK: nonlinear extension pair distances=' + str(dict(sorted(short_extension_histogram.items()))))
print('ADDITIVE BLOCKING: all 1395 dimension-3 subspaces meet the minimum-weight set')

def permute(x, p):
    return sum(((x >> i) & 1) << p[i] for i in range(6))

seed = next(iter(planes))
even_permutations = [p for p in permutations(range(6)) if sum(p[i] > p[j] for i in range(6) for j in range(i + 1, 6)) % 2 == 0]
three_cover = {frozenset(permute(x, p) for x in seed) for p in even_permutations}
assert len(three_cover) == 15
assert all(sum(x in w for w in three_cover) == (3 if x in F else 0) for x in range(1, 64))
coordinate_for_plane = {null_plane(m, e): (m, e) for m, e in source_coordinates}
three_coordinates = [coordinate_for_plane[w] for w in sorted(three_cover, key=lambda w: tuple(sorted(w)))]
for lam in range(2, 33):
    b = lam % 2
    a = (lam - 3 * b) // 2
    assert a >= 0 and 2 * a + 3 * b == lam
    coords = two_coordinates * a + three_coordinates * b
    assert len(coords) == 5 * lam
    assert all(sum(rho(m, e, x) == 0 for m, e in coords) == (lam if x in F else 0) for x in range(1, 64))
    words = code(coords)
    assert distances(words) == {4 * lam: 1120, 5 * lam: 896}
    assert min(distances(extension(words, source_tile_plane))) == 4 * lam + 1
print('FAMILY CHECK: 2-cover has 10 planes; 3-cover has 15 planes')
print('FAMILY CHECK: lambda=2..32; n=5*lambda; d=4*lambda; extension d+1; faithful')

# A 1-cover would require disjoint planes; exhaustive independent check finds
# no three pairwise disjoint planes, hence no partition into five planes.
disjoint_pairs = [(u, w) for u, w in combinations(planes, 2) if u & w == {0}]
assert disjoint_pairs
assert not any(u & z == {0} and w & z == {0} for u, w in disjoint_pairs for z in planes)
assert all(((x.bit_count() * (x.bit_count() + 1) // 2) % 2 == 0) == (x == 0 or x in F) for x in V)
assert all(sum((x.bit_count() * (x.bit_count() + 1) // 2) % 2 for x in w) % 2 == 0 for w in three_spaces)
print('CONSTRUCTION LIMIT: no three disjoint singular planes; no positive 1-cover')
print('SCOPE: finite checks support the outline; they do not prove the infinite family or global minimality of 10')
print('ALL_CHECKS_OK')
