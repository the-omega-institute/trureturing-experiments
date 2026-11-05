#!/usr/bin/env python3
"""Exact finite checks for Ivanov Problem 8.1; not an unbounded proof.

Requires Python 3 and NetworkX 3.6.1. All graph tests use integer labels,
integer degrees, and exact adjacency; no floating point is used as evidence.
"""
import itertools
import json
import networkx as nx


def clock():
    q = nx.Graph()
    q.add_nodes_from(range(12))
    for i in range(6):
        q.add_edge(i, (i + 1) % 6)
        for d in (0, 1, 3):
            q.add_edge(6 + i, (i + d) % 6)
    return q


def is_iso(g, h, f):
    if set(f) != set(g) or set(f.values()) != set(h) or len(f) != len(set(f.values())):
        return False
    return all(g.has_edge(u, v) == h.has_edge(f[u], f[v])
               for u, v in itertools.combinations(g, 2))


def auts(g):
    return list(nx.algorithms.isomorphism.GraphMatcher(g, g).isomorphisms_iter())


def without(g, v):
    return g.subgraph(set(g) - {v}).copy()


def assert_asymmetric(g):
    maps = auts(g)
    assert len(maps) == 1 and all(maps[0][v] == v for v in g)


def brute_core(q):
    # Degree classes force A and B to be invariant. Distinct B neighborhoods
    # force the B permutation once the A permutation is chosen. Enumerate all
    # 6! A permutations, without assuming a dihedral or rotational action.
    assert [q.degree(i) for i in range(12)] == [5] * 6 + [3] * 6
    k = nx.complement(q)
    assert [k.degree(i) for i in range(12)] == [6] * 6 + [8] * 6
    b_by_neighbors = {frozenset(q.neighbors(6 + i)): 6 + i for i in range(6)}
    assert len(b_by_neighbors) == 6
    candidates = []
    tested = 0
    for p in itertools.permutations(range(6)):
        tested += 1
        f = dict(enumerate(p))
        for i in range(6):
            image = frozenset(p[j] for j in q.neighbors(6 + i))
            if image not in b_by_neighbors:
                break
            f[6 + i] = b_by_neighbors[image]
        else:
            if is_iso(k, k, f):
                assert is_iso(q, q, f)
                candidates.append(tuple(f[i] for i in range(12)))
    rotations = {tuple((i + r) % 6 if i < 6 else 6 + (i - 6 + r) % 6
                       for i in range(12)) for r in range(6)}
    assert tested == 720 and set(candidates) == rotations
    assert {tuple(f[i] for i in range(12)) for f in auts(k)} == rotations
    assert nx.is_connected(k)
    return k, {"A_permutations_tested": tested, "automorphisms": len(candidates),
               "automorphism_maps": sorted(candidates), "K_connected": True}


def fidelity(q):
    # Source Section 3, Example 3.1: five leaves on Q_6 after x_0 is removed.
    # Its G_1 has 17 vertices, 29 edges, and the four specified active deletions.
    g = q.copy()
    g.add_edges_from((12 + i, i) for i in range(1, 6))
    assert len(g) == 17 and g.number_of_edges() == 29
    orders = []
    for i in range(1, 5):
        assert nx.is_connected(g)
        assert_asymmetric(g)
        a, b = 12 + i, 17
        assert a != b and a in g and b in g
        assert nx.is_isomorphic(without(g, a), without(g, b))
        orders.append(len(g))
        g.remove_node(a)
    assert nx.is_connected(g)
    assert_asymmetric(g)
    return {"source_example": "3.1", "starting_order": 17, "starting_edges": 29,
            "active_deletions": 4, "predeletion_orders": orders,
            "terminal_order": len(g)}


def family(k, m):
    g = k.copy()
    g.add_nodes_from(range(12, 12 + m))
    g.add_edges_from((12 + i, 12 + i + 1) for i in range(m - 1))
    g.add_edges_from((12 + i, i % 6) for i in range(m))
    return g


def check_family(k):
    data = []
    for m in range(1, 101):
        g = family(k, m)
        assert len(g) == 12 + m and nx.is_connected(g)
        assert max(g.degree(i) for i in range(12, 12 + m)) <= 3
        assert min(g.degree(i) for i in range(12)) >= 6
        assert_asymmetric(g)
        if m >= 2:
            left, right = 12, 12 + m - 1
            f = {i: (i + 1) % 6 for i in range(6)}
            f.update({6 + i: 6 + (i + 1) % 6 for i in range(6)})
            f.update({12 + i: 12 + i + 1 for i in range(m - 1)})
            assert is_iso(without(g, right), without(g, left), f)
            assert nx.is_isomorphic(without(g, right), without(g, left))
        data.append({"m": m, "order": len(g), "connected": True,
                     "automorphisms": 1, "endpoint_cards_isomorphic": m >= 2})
    return data


def check_all_stages(k):
    results = []
    for n in range(14, 19):
        m = n - 12
        g = family(k, m)
        stages = []
        for remaining in range(m, 0, -1):
            assert nx.is_connected(g)
            assert_asymmetric(g)
            assert set(g) == set(family(k, remaining))
            assert set(map(frozenset, g.edges())) == set(map(frozenset, family(k, remaining).edges()))
            stage = {"order": len(g), "connected": True, "automorphisms": 1}
            if remaining >= 2:
                deleted, mate = 12 + remaining - 1, 12
                assert deleted != mate
                assert nx.is_isomorphic(without(g, deleted), without(g, mate))
                assert all(f[deleted] != mate for f in auts(g))
                stage.update({"deleted": deleted, "mate": mate, "pseudosimilar": True})
                g.remove_node(deleted)
            stages.append(stage)
        assert len(stages) - 1 == n - 13
        results.append({"starting_order": n, "chain_length": n - 13, "stages": stages})
    return results


def main():
    q = clock()
    k, core = brute_core(q)
    result = {"networkx_version": nx.__version__, "fidelity": fidelity(q),
              "core": core, "family_m_1_through_100": check_family(k),
              "all_stages_n_14_through_18": check_all_stages(k),
              "scope": "Finite corroboration only; does not prove the unbounded statement."}
    print(json.dumps(result, indent=2, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
