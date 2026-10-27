#!/usr/bin/env python3
"""Exact finite checks for Schachner--Udell, arXiv:2610.18463v1, Q5.9.

Permutations are integer tuples of images; composition p*q means p after q.
No symbolic family lemma is assumed by the enumerator. Python standard library only.
"""
from itertools import combinations, permutations
from math import factorial


def edge(u, v):
    return (min(u, v), max(u, v))


def hops(n, edges):
    """Enumerate all adjacency-respecting bijections by independent backtracking."""
    neighbors = [sorted(v for v in range(n) if edge(u, v) in edges) for u in range(n)]
    order = sorted(range(n), key=lambda u: len(neighbors[u]))
    image = [-1] * n
    found = set()

    def visit(k, used):
        if k == n:
            found.add(tuple(image))
            return
        u = order[k]
        for v in neighbors[u]:
            if not (used >> v) & 1:
                image[u] = v
                visit(k + 1, used | (1 << v))
        image[u] = -1

    visit(0, 0)
    return found


def brute_hops(n, edges):
    return {p for p in permutations(range(n)) if all(edge(u, p[u]) in edges for u in range(n))}


def compose(p, q):
    return tuple(p[q[u]] for u in range(len(p)))


def generated(n, generators):
    """Breadth-first closure by right multiplication, including the identity."""
    identity = tuple(range(n))
    group, queue = {identity}, [identity]
    for p in queue:
        for q in generators:
            r = compose(p, q)
            if r not in group:
                group.add(r)
                queue.append(r)
    # Independent verification of the computed subgroup, including inverses.
    for p in group:
        inverse = tuple(p.index(u) for u in range(n))
        assert inverse in group
    assert all(compose(p, q) in group for p in group for q in group)
    return group


def connected(n, edges):
    reached, queue = {0}, [0]
    for u in queue:
        for v in range(n):
            if edge(u, v) in edges and v not in reached:
                reached.add(v)
                queue.append(v)
    return len(reached) == n


def maximal(ambient, chosen, invariant):
    old = invariant(chosen)
    return all(old < invariant(chosen | {e}) for e in ambient - chosen)


def family(t):
    n = 4 + 2 * t
    ambient = set(combinations(range(4), 2))
    for i in range(t):
        u, v = 4 + 2 * i, 5 + 2 * i
        ambient |= {edge(0, u), edge(u, v)}
    h1, h2 = ambient - {edge(1, 3)}, ambient - {edge(0, 2)}
    return n, ambient, h1, h2, h1 & h2


def extend(core, t):
    return tuple(core) + tuple(v for i in range(t) for v in (5 + 2 * i, 4 + 2 * i))


def preserves_partition(p):
    return {p[0], p[2]} in ({0, 2}, {1, 3})


def main():
    # Source fidelity: Examples 1.2(i), Proposition 2.2, and the C8 computation.
    for n in range(3, 9):
        edges = {edge(u, (u + 1) % n) for u in range(n)}
        hs = hops(n, edges)
        assert hs == brute_hops(n, edges)
        assert len(hs) == (2 if n % 2 else 4)
        group = generated(n, hs)
        assert len(group) == (n if n % 2 else n * n // 2)
        print(f'FIDELITY C{n}: hops={len(hs)} leap_order={len(group)}')
    for n in range(2, 7):
        edges = set(combinations(range(n), 2))
        hs = hops(n, edges)
        assert hs == brute_hops(n, edges)
        group = generated(n, hs)
        assert len(group) == (3 if n == 3 else factorial(n))
        print(f'FIDELITY K{n}: hops={len(hs)} leap_order={len(group)}')

    core_hops = {(1, 2, 3, 0), (3, 0, 1, 2), (1, 0, 3, 2), (3, 2, 1, 0)}
    for t in range(7):
        n, ambient, h1, h2, common = family(t)
        graphs = (common, h1, h2, ambient)
        hs = [hops(n, e) for e in graphs]
        if t <= 2:
            assert all(h == brute_hops(n, e) for h, e in zip(hs, graphs))
        groups = [generated(n, h) for h in hs]
        assert [len(e) for e in graphs] == [4 + 2*t, 5 + 2*t, 5 + 2*t, 6 + 2*t]
        assert all(connected(n, e) for e in graphs)
        expected = {extend(p, t) for p in core_hops}
        assert hs[0] == hs[1] == hs[2] == expected
        assert [len(h) for h in hs] == [4, 4, 4, 9]
        assert groups[0] == groups[1] == groups[2]
        assert [len(g) for g in groups] == [8, 8, 8, 24 if t == 0 else 48]
        assert all(preserves_partition(p) for p in groups[0])
        witness = extend((2, 3, 1, 0), t)  # (1 3 2 4) tau in one-based labels
        assert witness in hs[3] and witness not in groups[0]
        assert not preserves_partition(witness)
        for h in hs[3]:
            assert all(h[4+2*i] == 5+2*i and h[5+2*i] == 4+2*i for i in range(t))
        hop_invariant = lambda e: hops(n, e)
        leap_invariant = lambda e: generated(n, hops(n, e))
        for invariant in (hop_invariant, leap_invariant):
            assert maximal(ambient, h1, invariant)
            assert maximal(ambient, h2, invariant)
            assert not maximal(ambient, common, invariant)
        print(f'FAMILY t={t} n={n}: connected=4/4 edges={[len(e) for e in graphs]} '
              f'hops={[len(h) for h in hs]} leap_orders={[len(g) for g in groups]} '
              'both_intersection_obstructions=True')

    # All spanning subgraphs of K4: exhibit failure for the complete families.
    n, ambient, h1, h2, common = family(0)
    edges = sorted(ambient)
    subsets = [{edges[i] for i in range(6) if (mask >> i) & 1} for mask in range(64)]
    for label, invariant in [('hop', lambda e: hops(4, e)),
                              ('leap', lambda e: generated(4, hops(4, e)))]:
        family_of_sets = {frozenset(e) for e in subsets if maximal(ambient, e, invariant)}
        assert frozenset(h1) in family_of_sets and frozenset(h2) in family_of_sets
        assert frozenset(common) not in family_of_sets
        failures = [(x, y) for x in family_of_sets for y in family_of_sets if x & y not in family_of_sets]
        assert failures
        print(f'K4 ALL_SUBGRAPHS {label}: checked=64 maximal_sets={len(family_of_sets)} '
              f'ordered_intersection_failures={len(failures)}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
