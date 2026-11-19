"""Exact checks for Aboulker et al., arXiv:2310.04265v2, Question 3.5."""

from itertools import combinations, permutations, product


def edges(pairs):
    return {tuple(sorted(pair)) for pair in pairs}


def backedges(arc, order):
    return edges((earlier, later) for earlier, later in combinations(order, 2)
                 if arc(later, earlier))


def adjacency(vertices, graph):
    neighbors = {vertex: set() for vertex in vertices}
    for u, v in graph:
        neighbors[u].add(v)
        neighbors[v].add(u)
    return neighbors


def triangle_free(graph, neighbors):
    return all(not (neighbors[u] & neighbors[v]) for u, v in graph)


def colorable(neighbors, limit):
    """Exhaustive DSATUR search, trying unused colors up to relabeling."""
    colors = {}
    uncolored = set(neighbors)

    def search(used):
        if not uncolored:
            return True
        forbidden = {
            v: {colors[u] for u in neighbors[v] if u in colors}
            for v in uncolored
        }
        v = max(uncolored, key=lambda v: (len(forbidden[v]), len(neighbors[v]), -v))
        uncolored.remove(v)
        for c in range(min(used + 1, limit)):
            if c not in forbidden[v]:
                colors[v] = c
                if search(max(used, c + 1)):
                    return True
                del colors[v]
        uncolored.add(v)
        return False

    return search(0)


def is_forward(arc, order):
    return all(arc(u, v) for u, v in combinations(order, 2))


def check_fidelity():
    """The example after Theorem 3.3 uses increasing and reversed orders of TT_n."""
    for n in range(1, 7):
        vertices = range(n)
        readings = {}
        for order in permutations(vertices):
            graph = backedges(int.__lt__, order)
            clique = max((len(group) for r in range(n + 1)
                          for group in combinations(vertices, r)
                          if edges(combinations(group, 2)) <= graph))
            neighbors = adjacency(vertices, graph)
            chromatic = next(k for k in range(1, n + 1) if colorable(neighbors, k))
            readings[order] = (clique, chromatic)
        assert readings[tuple(reversed(vertices))] == (n, n)
        assert tuple(map(min, zip(*readings.values()))) == (1, 1)
    print("Fidelity TT_n, n=1..6: reverse backedge (omega, chi)=(n,n); minima=(1,1)")


def separator():
    a, b, c, x, y, z, u = range(7)
    left, right = (a, b, c), (x, y, z)
    arcs = {(a, b), (b, c), (c, a), (x, y), (y, z), (z, x)}
    arcs.update(product(left, right))
    arcs.update((u, v) for v in left)
    arcs.update((v, u) for v in right)

    def arc(s, t):
        return (s, t) in arcs

    order = (a, b, x, u, c, y, z)
    color_orders = ((u,), (a, b, x, y), (c, z))
    assert all(arc(s, t) != arc(t, s) for s, t in combinations(range(7), 2))
    assert not any(arc(v, v) for v in range(7))
    assert all(is_forward(arc, part) for part in color_orders)
    graph = backedges(arc, order)
    assert graph == edges(((u, a), (u, b), (u, y), (u, z), (a, c), (c, x), (x, z)))
    assert triangle_free(graph, adjacency(range(7), graph))

    cyclic_triples = [
        (s, t, v) for s, t, v in combinations(range(7), 3)
        if arc(s, t) == arc(t, v) == arc(v, s)
    ]
    assert all(any(colors[s] == colors[t] == colors[v] for s, t, v in cyclic_triples)
               for colors in product(range(2), repeat=7))

    # The endpoints are 7 and 8. In every tested order, endpoint 8 precedes endpoint 7.
    # All arcs through the separator have the direction 7 -> S -> 8.
    # The conclusion must hold for either arc between endpoints.
    for forward_endpoint_arc in (False, True):
        extended_arcs = (
            arcs
            | {(7, v) for v in range(7)}
            | {(v, 8) for v in range(7)}
            | ({(7, 8)} if forward_endpoint_arc else {(8, 7)})
        )
        tested = 0
        for extended_order in permutations(range(9)):
            if extended_order.index(7) < extended_order.index(8):
                continue
            graph = backedges(lambda s, t: (s, t) in extended_arcs, extended_order)
            assert not triangle_free(graph, adjacency(range(9), graph))
            tested += 1
        assert tested == 181440
    print("Separator: 7 edges in backedge graph, d=3; all 362880 reversed-endpoint orders obstructed")
    return arc, order, color_orders


def mycielski(n, graph):
    doubled = edges((u, v + n) for a, b in graph for u, v in ((a, b), (b, a)))
    apex_edges = {(v, 2 * n) for v in range(n, 2 * n)}
    return 2 * n + 1, graph | doubled | apex_edges


def check_family(q, n, graph, colors_without_apex, gadget):
    gadget_arc, gadget_order, gadget_color_orders = gadget
    assert n == 6 * 2 ** (q - 3) - 1
    assert all(v == n - 1 or colors_without_apex[u] != colors_without_apex[v]
               for u, v in graph)
    neighbors = adjacency(range(n), graph)
    assert graph and triangle_free(graph, neighbors)
    assert not colorable(neighbors, q - 1)
    assert colorable(neighbors, q)

    apex = n - 1
    first_class = [v for v, c in enumerate(colors_without_apex) if c == 0]
    others = [v for v, c in enumerate(colors_without_apex) if c != 0]
    base_order = ([v for v in first_class if v not in neighbors[apex]]
                  + [apex] + [v for v in first_class if v in neighbors[apex]] + others)
    base_rank = {v: rank for rank, v in enumerate(base_order)}

    def piece(v):
        return 2 * base_rank[v] if v < n else 2 * ((v - n) // 7) + 1

    def arc(u, v):
        if u < n and v < n:
            forward = base_rank[u] < base_rank[v]
            return forward != (tuple(sorted((u, v))) in graph)
        if piece(u) == piece(v):
            return gadget_arc((u - n) % 7, (v - n) % 7)
        return piece(u) < piece(v)

    total = n + 7 * (n - 1)
    assert total == 6 * 2 ** q - 15
    order = []
    expected = set(graph)
    gadget_graph = backedges(gadget_arc, gadget_order)
    for i, v in enumerate(base_order):
        order.append(v)
        if i < n - 1:
            start = n + 7 * i
            order.extend(start + s for s in gadget_order)
            expected.update((start + s, start + t) for s, t in gadget_graph)
            assert all(arc(v, start + s) and arc(start + s, base_order[i + 1])
                       for s in range(7))

    assert sorted(order) == list(range(total))
    assert all(arc(u, v) != arc(v, u) for u, v in combinations(order, 2))
    assert not any(arc(v, v) for v in order)
    observed = backedges(arc, order)
    assert observed == expected
    assert triangle_free(observed, adjacency(range(total), observed))

    color_orders = [[] for _ in range(q - 1)]
    for v in order:
        if v < n and v != apex:
            color_orders[colors_without_apex[v]].append(v)
        elif v >= n:
            local = (v - n) % 7
            color = 1 + next(i for i, group in enumerate(gadget_color_orders)
                             if local in group)
            color_orders[color].append(v)
    color_orders[0].append(apex)
    assert sorted(v for group in color_orders for v in group) == list(range(total))
    assert all(is_forward(arc, group) for group in color_orders)
    print(f"q={q}: base={n}, tournament={total}, directed omega=2, "
          f"d<={q - 1}, min chi over clique-optimal orders={q}")


def main():
    if not __debug__:
        raise RuntimeError("Run without -O so that assertions execute")
    check_fidelity()
    gadget = separator()
    n = 5
    graph = edges((v, (v + 1) % n) for v in range(n))
    colors = [0, 1, 0, 1, 2]
    for q in range(4, 7):
        n, graph = mycielski(n, graph)
        colors_without_apex = colors + colors
        if q >= 5:
            check_family(q, n, graph, colors_without_apex, gadget)
        colors = colors_without_apex + [q - 1]
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
