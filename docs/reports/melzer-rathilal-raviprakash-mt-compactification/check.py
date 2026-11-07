#!/usr/bin/env python3
"""Exact finite checks for arXiv:2606.05832v1, Remark 3.10.

An integer mask represents a subset. All topologies on labelled sets of
zero through three points, and all frame homomorphisms between their
open-set frames, are included. This is finite evidence, not an unbounded
proof or a test of representation of nonspatial frames.
"""

from itertools import product


def join(items):
    result = 0
    for item in items:
        result |= item
    return result


def le(a, b):
    return a & b == a


def topologies(n):
    top = (1 << n) - 1
    optional = [a for a in range(1 << n) if a not in (0, top)]
    for selection in range(1 << len(optional)):
        opens = tuple(sorted({0, top} | {
            a for i, a in enumerate(optional) if selection & (1 << i)
        }))
        members = set(opens)
        if all(a & b in members and a | b in members for a, b in product(opens, repeat=2)):
            # In a finite set, binary closure and the empty union suffice.
            assert all(join(a for i, a in enumerate(opens) if selection & (1 << i))
                       in members for selection in range(1 << len(opens)))
            yield opens


def interior(opens, carrier_size):
    return tuple(join(u for u in opens if le(u, x)) for x in range(carrier_size))


def is_frame_hom(source, target, values):
    g = dict(zip(source, values))
    return (g[0] == 0 and g[source[-1]] == target[-1]
            and all(g[a & b] == g[a] & g[b] and g[a | b] == g[a] | g[b]
                    for a, b in product(source, repeat=2)))


def frame_homs(source, target):
    irreducibles = [a for a in source if a != 0 and all(
        b | c != a or b == a or c == a for b, c in product(source, repeat=2))]
    below = {a: [i for i, b in enumerate(irreducibles) if le(b, a)] for a in source}
    assert all(join(irreducibles[i] for i in below[a]) == a for a in source)
    for images in product(target, repeat=len(irreducibles)):
        values = tuple(join(images[i] for i in below[a]) for a in source)
        if is_frame_hom(source, target, values):
            yield values


def verify_projection(n, m):
    size = 1 << (n + m)
    ctop, mtop = size - 1, (1 << m) - 1
    p = lambda x: x & mtop
    assert p(0) == 0 and p(ctop) == mtop
    assert set(map(p, range(size))) == set(range(1 << m))
    for a in range(size):
        assert p(ctop ^ a) == mtop ^ p(a)
        for b in range(size):
            assert p(a & b) == p(a) & p(b)
            assert p(a | b) == p(a) | p(b)


def fidelity():
    # Source Example 3.8: a two-element discrete Boolean algebra maps
    # into the four-element indiscrete one. Its opens map is an
    # isomorphism, but it is not an MT-compactification: not onto.
    nbox = (0, 1)
    mbox = interior((0, 3), 4)
    h = (0, 3)
    assert all(le(h[nbox[a]], mbox[h[a]]) for a in range(2))
    assert all(h[1 ^ a] == 3 ^ h[a] for a in range(2))
    assert all(h[a & b] == h[a] & h[b] and h[a | b] == h[a] | h[b]
               for a, b in product(range(2), repeat=2))
    assert {h[a] for a in range(2) if nbox[a] == a} == {0, 3}
    assert set(h) != set(range(4))
    assert is_frame_hom((0, 1), (0, 3), h)
    print('FIDELITY Example 3.8: opens isomorphism, Boolean carrier not onto: OK')


def check_realization(n, m, source, target, values):
    g = dict(zip(source, values))
    size, mtop = 1 << (n + m), (1 << m) - 1
    ctop = size - 1
    graph = {a: g[a] | (a << m) for a in source}
    graph_set = set(graph.values())
    box = interior(graph_set, size)
    mbox = interior(target, 1 << m)
    assert len(graph_set) == len(source)
    assert graph[0] == 0 and graph[source[-1]] == ctop
    assert box[ctop] == ctop
    assert {a for a in range(size) if box[a] == a} == graph_set
    for a in range(size):
        assert le(box[a], a) and box[box[a]] == box[a]
        assert box[a] in graph_set
        assert all(le(s, box[a]) for s in graph_set if le(s, a))
        assert le(box[a] & mtop, mbox[a & mtop])
        for b in range(size):
            assert box[a & b] == box[a] & box[b]
    for a, b in product(source, repeat=2):
        assert le(graph[a], graph[b]) == le(a, b)
        assert graph[a & b] == graph[a] & graph[b]
        assert graph[a | b] == graph[a] | graph[b]
    assert all(graph[a] & mtop == g[a] for a in source)
    # Enumerate every cover by open elements; each cover is its own
    # finite subcover, including the empty cover in the trivial case.
    for selection in range(1 << len(source)):
        subset = [a for i, a in enumerate(source) if selection & (1 << i)]
        assert graph[join(subset)] == join(graph[a] for a in subset)
        if join(graph[a] for a in subset) == ctop:
            assert join(subset) == source[-1]
    dense = all(g[a] != 0 or a == 0 for a in source)
    onto = set(values) == set(target)
    if dense and onto:
        assert all(s & mtop != 0 or s == 0 for s in graph_set)
        assert {s & mtop for s in graph_set} == set(target)
        assert {a & mtop for a in range(size)} == set(range(1 << m))
    return int(dense and onto)


def main():
    fidelity()
    all_topologies = {n: tuple(topologies(n)) for n in range(4)}
    assert [len(all_topologies[n]) for n in range(4)] == [1, 1, 4, 29]
    print('TOPOLOGIES labelled point counts 0,1,2,3: 1,1,4,29')
    for n, families in all_topologies.items():
        for opens in families:
            box = interior(opens, 1 << n)
            assert {a for a in range(1 << n) if box[a] == a} == set(opens)
    # Independently enumerate all function tables for the smaller cases.
    smaller_pairs = 0
    for n, m in product(range(3), repeat=2):
        for source, target in product(all_topologies[n], all_topologies[m]):
            brute = {values for values in product(target, repeat=len(source))
                     if is_frame_hom(source, target, values)}
            assert set(frame_homs(source, target)) == brute
            smaller_pairs += 1
    print('ENUMERATOR independent full-function cross-check frame pairs:', smaller_pairs)
    homs = compacts = pairs = 0
    for n, m in product(range(4), repeat=2):
        verify_projection(n, m)
        count = compact_count = 0
        for source, target in product(all_topologies[n], all_topologies[m]):
            pairs += 1
            for values in frame_homs(source, target):
                compact_count += check_realization(n, m, source, target, values)
                count += 1
        homs += count
        compacts += compact_count
        print('GRAPH point counts', n, m, 'frame_homs', count,
              'dense_onto_compactifications', compact_count)
    print('TOTAL frame_pairs', pairs, 'frame_homs', homs, 'compactifications', compacts)
    print('LIMIT finite spatial representatives only; no proof for arbitrary frames')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
