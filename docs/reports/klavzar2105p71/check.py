#!/usr/bin/env python3
"""Exact checks of Steiner general position for finite simple graphs."""

from itertools import combinations


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def mask(vertices):
    return sum(1 << vertex for vertex in set(vertices))


def members(vertex_mask):
    while vertex_mask:
        low = vertex_mask & -vertex_mask
        yield low.bit_length() - 1
        vertex_mask ^= low


def adjacency(order, edges):
    neighbors = [0] * order
    for u, v in edges:
        require(0 <= u < order and 0 <= v < order and u != v, (u, v))
        neighbors[u] |= 1 << v
        neighbors[v] |= 1 << u
    return neighbors


def connected(neighbors, allowed):
    if not allowed:
        return False
    reached = frontier = allowed & -allowed
    while frontier:
        candidates = 0
        for vertex in members(frontier):
            candidates |= neighbors[vertex]
        frontier = candidates & allowed & ~reached
        reached |= frontier
    return reached == allowed


def steiner_supports(neighbors, k):
    """For each terminal set, store minimum edges and union of minimum tree supports.

    A vertex set is the support of a tree exactly when its induced graph
    is connected. Every tree on it has cardinality minus one edges.
    """
    order = len(neighbors)
    records = {}
    for size in range(k, order + 1):
        for vertices in combinations(range(order), size):
            support = mask(vertices)
            if not connected(neighbors, support):
                continue
            for terminals in combinations(vertices, k):
                terminal_mask = mask(terminals)
                best_size, best_support = records.get(terminal_mask, (size, 0))
                if size == best_size:
                    records[terminal_mask] = (size, best_support | support)
        # The entire cardinality layer is needed to capture ties.
        if all(mask(terms) in records for terms in combinations(range(order), k)):
            return {
                terminals: (size - 1, support)
                for terminals, (size, support) in records.items()
            }
    raise ValueError("expected a connected host graph and 1 <= k <= order")


def forbidden_sets(records):
    return {
        terminals | (1 << extra)
        for terminals, (_, support) in records.items()
        for extra in members(support & ~terminals)
    }


def feasible(forbidden, selected):
    return all(selected & obstruction != obstruction for obstruction in forbidden)


def maximum_position(order, forbidden):
    for size in range(order, -1, -1):
        for vertices in combinations(range(order), size):
            if feasible(forbidden, mask(vertices)):
                return size
    raise AssertionError("the empty set is feasible")


def position_number(neighbors, k):
    return maximum_position(len(neighbors), forbidden_sets(steiner_supports(neighbors, k)))


def check_source():
    cases = 0
    for order in range(3, 10):
        path_edges = [(v, v + 1) for v in range(order - 1)]
        for edges in (path_edges, [(0, v) for v in range(1, order)]):
            tree = adjacency(order, edges)
            leaves = sum(degree_mask.bit_count() == 1 for degree_mask in tree)
            for k in range(2, order):
                actual = position_number(tree, k)
                require(actual == max(leaves, k), ("Theorem 3.1", order, leaves, k, actual))
                cases += 1
        cycle = adjacency(order, [*path_edges, (order - 1, 0)])
        for k in range(2, order):
            expected = k if 2 * order // 3 <= k <= order - 2 else k + 1
            actual = position_number(cycle, k)
            require(actual == expected, ("Theorem 3.2", order, k, actual, expected))
            cases += 1
    print(f"Fidelity: Theorems 3.1 and 3.2, {cases} exact path/star/cycle cases")

    join = adjacency(7, [(0, 1), *((u, v) for u in (0, 1) for v in range(2, 7))])
    actual = tuple(position_number(join, k) for k in (2, 3))
    require(actual == (5, 5), ("complete split graph", actual))
    print(f"Corollary 4.4(v) discrepancy: K_2 join independent 5, k=2,3 -> {actual}")


def check_family():
    names = {letter: index for index, letter in enumerate("ABCDEFpq")}
    certificate_rows = (
        ("ABCD", "D", "AC CD DB"), ("ACEF", "A", "AC AE AF"),
        ("BDEF", "B", "BD BE BF"), ("ABCp", "A", "pA AC AF FB"),
        ("ACEp", "A", "pA AC AE"), ("ADFp", "A", "pA AF FB BD"),
        ("BDEp", "E", "pE EB BD"), ("BEFp", "E", "pE EB EF"),
        ("CDFp", "C", "pA AC CD AF"), ("ACpq", "A", "CA Ap Aq"),
        ("ADpq", "A", "DC CA Ap Aq"), ("AFpq", "A", "FA Ap Aq"),
        ("BDpq", "B", "DB BE Ep Eq"), ("BEpq", "E", "BE Ep Eq"),
        ("CDpq", "C", "DC CA Ap Aq"), ("DEpq", "E", "DB BE Ep Eq"),
        ("EFpq", "E", "FE Ep Eq"), ("BCFp", "F", "pA AC AF FB"),
    )
    for m in range(2, 9):
        order = m + 6
        core = [tuple(names[x] for x in edge) for edge in "AC BD CD AE BE AF BF EF".split()]
        edges = [*core, *((u, v) for u in (0, 4) for v in range(6, order))]
        graph = adjacency(order, edges)
        all_vertices = (1 << order) - 1
        require(connected(graph, all_vertices), ("connected", m))
        require(
            all(connected(graph, all_vertices ^ (1 << v)) for v in range(order)),
            ("2-connected", m),
        )
        pair_records = steiner_supports(graph, 2)
        triple_records = steiner_supports(graph, 3)
        pair_obstructions = forbidden_sets(pair_records)
        triple_obstructions = forbidden_sets(triple_records)
        twins = list(range(6, order))
        two_witness = mask([*twins, *(names[x] for x in "BCF")])
        three_witness = mask([*twins, *(names[x] for x in "CF")])
        require(feasible(pair_obstructions, two_witness), ("pair lower witness", m))
        require(feasible(triple_obstructions, three_witness), ("triple lower witness", m))
        actual = tuple(maximum_position(order, obs) for obs in (pair_obstructions, triple_obstructions))
        require(actual == (m + 3, m + 2), ("family maxima", m, actual))
        for four, extra, tree in certificate_rows:
            selected = mask(names[x] for x in four)
            terminal = selected ^ (1 << names[extra])
            tree_edges = [tuple(names[x] for x in edge) for edge in tree.split()]
            support = mask(x for edge in tree_edges for x in edge)
            require(selected & support == selected, ("certificate support", m, four))
            require(len(tree_edges) == support.bit_count() - 1, ("certificate size", m, four))
            require(connected(adjacency(order, tree_edges), support), ("certificate tree", m, four))
            require(all(graph[u] & (1 << v) for u, v in tree_edges), ("certificate edges", m, four))
            require(triple_records[terminal][0] == len(tree_edges), ("certificate optimum", m, four))
        late = position_number(graph, order - 1)
        require(late > actual[1], ("later increase", m, late))
        print(f"G_{m}: order={order}, sgp_2={actual[0]}, sgp_3={actual[1]}, sgp_{order - 1}={late}, 18 certificates")


if __name__ == "__main__":
    check_source()
    check_family()
    print("ALL_CHECKS_OK")
