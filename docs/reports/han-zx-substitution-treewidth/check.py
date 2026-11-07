#!/usr/bin/env python3
"""Exact checks for Han, arXiv:2611.01890v1, Conjecture A.8.

Standard library only. Elimination uses integer bitsets and a full memoized
minimum over choices; the independent bag checker uses the source definition.
"""
from functools import lru_cache
from itertools import combinations
import json


def graph(vertices, edges):
    vertices = frozenset(vertices)
    edges = frozenset(tuple(sorted(e)) for e in edges)
    assert all(u in vertices and v in vertices and u != v for u, v in edges)
    return vertices, edges


def clique_edges(vertices):
    return set(combinations(sorted(vertices), 2))


def eliminate(adj, alive, v):
    neighbors = adj[v] & alive
    rest = alive & ~(1 << v)
    new = list(adj)
    for u in range(len(adj)):
        if neighbors & (1 << u):
            new[u] = (new[u] | neighbors) & rest & ~(1 << u)
    new[v] = 0
    return tuple(new), rest, neighbors.bit_count()


def exact_treewidth(g):
    vertices, edges = g
    labels = sorted(vertices)
    index = {v: i for i, v in enumerate(labels)}
    adj = [0] * len(labels)
    for u, v in edges:
        i, j = index[u], index[v]
        adj[i] |= 1 << j
        adj[j] |= 1 << i

    @lru_cache(None)
    def best(adj, alive):
        if not alive:
            return -1
        return min(max(degree, best(new, rest))
                   for v in range(len(adj)) if alive & (1 << v)
                   for new, rest, degree in [eliminate(adj, alive, v)])

    adj, alive = tuple(adj), (1 << len(labels)) - 1
    width = best(adj, alive)
    order = []
    while alive:
        target = best(adj, alive)
        for v in range(len(labels)):
            if not alive & (1 << v):
                continue
            new, rest, degree = eliminate(adj, alive, v)
            if max(degree, best(new, rest)) == target:
                order.append(labels[v])
                adj, alive = new, rest
                break
        else:
            raise AssertionError('No minimizing elimination choice')
    return width, order, best.cache_info().currsize


def check_elimination(g, order):
    vertices, edges = g
    assert set(order) == vertices and len(order) == len(vertices)
    neighbors = {v: set() for v in vertices}
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    width = -1
    for v in order:
        ns = neighbors[v]
        width = max(width, len(ns))
        for u in ns:
            neighbors[u].update(ns - {u})
            neighbors[u].discard(v)
        del neighbors[v]
    return width


def check_bags(g, bags):
    vertices, edges = g
    bags = [set(b) for b in bags]
    assert bags and set().union(*bags) == vertices
    assert all(b <= vertices for b in bags)
    assert all(any(u in b and v in b for b in bags) for u, v in edges)
    for v in vertices:
        indices = [i for i, b in enumerate(bags) if v in b]
        assert indices == list(range(min(indices), max(indices) + 1))
    return max(map(len, bags)) - 1


def connected(g, deleted=None):
    vertices, edges = g
    vertices = vertices - {deleted}
    if not vertices:
        return False
    seen = {min(vertices)}
    while True:
        new = seen | {v for u, v in edges if u in seen and v in vertices}
        new |= {u for u, v in edges if v in seen and u in vertices}
        if new == seen:
            return seen == vertices
        seen = new


def complete_loose(core, endpoints, leaves):
    vertices, edges = core
    assert len(endpoints) == len(leaves) == len(set(leaves))
    assert not vertices.intersection(leaves)
    assert set(endpoints) <= vertices
    return graph(vertices | set(leaves), set(edges) | set(zip(endpoints, leaves)))


def substitute(g, region, replacement, ports):
    vertices, edges = g
    new_vertices, new_edges = replacement
    outside = vertices - region
    cut = {(u, v) if u in region else (v, u)
           for u, v in edges if (u in region) != (v in region)}
    assert len(ports) == len(cut)
    assert {(a, b) for a, b, s in ports} == cut
    assert all(s in new_vertices for a, b, s in ports)
    assert not vertices.intersection(new_vertices)
    result = graph(outside | new_vertices,
                   {e for e in edges if set(e) <= outside} | set(new_edges) |
                   {(b, s) for a, b, s in ports})
    return result, cut


def family(r):
    assert r >= 4
    a, b = set(range(r)), set(range(r, 2 * r))
    s1, s2 = 2 * r, 2 * r + 1
    g = graph(a | b, clique_edges(a) | clique_edges(b) | {(0, r), (1, r + 1)})
    h_core = graph(a, {e for e in g[1] if set(e) <= a})
    s_core = graph({s1, s2}, {(s1, s2)})
    gp, cut = substitute(g, a, s_core, [(0, r, s1), (1, r + 1, s2)])
    h = complete_loose(h_core, [0, 1], [2 * r, 2 * r + 1])
    s = complete_loose(s_core, [s1, s2], [2 * r + 2, 2 * r + 3])
    gs = [g, h, s, gp]
    bags = [[sorted(a), [0, 1, r, r + 1], sorted(b)],
            [[0, 2 * r], sorted(a), [1, 2 * r + 1]],
            [[2 * r + 2, s1], [s1, s2], [s2, 2 * r + 3]],
            [sorted(b), [r, r + 1, s1], [r + 1, s1, s2]]]
    widths = [check_bags(g, d) for g, d in zip(gs, bags)]
    assert widths == [r - 1, r - 1, 1, r - 1]
    # Explicit clique certificates; their lower-bound meaning is the classical
    # clique-in-a-bag theorem, not a new kernel-verified assertion here.
    for g, clique in zip(gs, [a, a, {s1, s2}, b]):
        assert clique <= g[0] and clique_edges(clique) <= g[1]
    assert connected(h_core) and connected(s_core)
    assert len(cut) == 2
    assert len({x for x, y in cut}) == len({y for x, y in cut}) == 2
    for g in [g, gp]:
        assert len(g[0]) >= 3 and connected(g)
        assert all(connected(g, v) for v in g[0])
    return gs, bags, widths


def main():
    # Fidelity to source Example 3.8: complete graphs have treewidth n - 1.
    for n in range(1, 8):
        g = graph(range(n), combinations(range(n), 2))
        tw, order, states = exact_treewidth(g)
        assert tw == n - 1 == check_elimination(g, order)
        print(f'FIDELITY complete_graph n={n} exact_treewidth={tw} states={states}')
    # Fidelity to Example 3.5: every nontrivial tree has treewidth 1.
    for n in range(2, 9):
        for name, edges in [('path', [(i, i + 1) for i in range(n - 1)]),
                            ('star', [(0, i) for i in range(1, n)])]:
            g = graph(range(n), edges)
            tw, order, states = exact_treewidth(g)
            assert tw == 1 == check_elimination(g, order)
            print(f'FIDELITY {name} n={n} exact_treewidth={tw} states={states}')
    for r in range(4, 8):
        gs, bags, widths = family(r)
        exact = []
        orders = []
        for g in gs:
            tw, order, states = exact_treewidth(g)
            assert tw == check_elimination(g, order)
            exact.append(tw)
            orders.append(order)
        assert exact == widths
        rhs = exact[0] - exact[1] + exact[2] + 2
        assert rhs == 3
        assert (exact[3] > rhs) == (r >= 5)
        print(f'FAMILY r={r} exact_widths={exact} k=2 rhs={rhs} '
              f'lhs={exact[3]} violation={exact[3] - rhs} 2_connected=True')
        if r == 5:
            names = ['G', 'H_completed', 'S_completed', 'G_prime']
            instance = {name: {'vertices': sorted(g[0]), 'edges': sorted(g[1]),
                              'path_bags': bag, 'width': tw, 'elimination_order': order}
                        for name, g, bag, tw, order in zip(names, gs, bags, exact, orders)}
            print('INSTANCE_R5 ' + json.dumps(instance, sort_keys=True))
    for r in range(8, 41):
        _, _, widths = family(r)
        assert widths[3] - (widths[0] - widths[1] + widths[2] + 2) == r - 4
    print('STRUCTURAL family r=8..40: explicit_bags_cliques_boundary_and_2_connectivity_OK')
    print('UNBOUNDED_FAMILY: the argument for all r>=5 is in README.md; finite checks do not prove it')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
