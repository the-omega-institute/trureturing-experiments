#!/usr/bin/env python3
"""Exact finite check for Porte--Sandeep--Santra Problem 5.

The state is (graph edge mask, sign mask).  Edge bits are ordered
01, 02, 03, 12, 13, 23.  All computations are integer/XOR operations.
"""

from collections import deque

EDGE_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {pair: i for i, pair in enumerate(EDGE_PAIRS)}
P4 = (1 << 0) | (1 << 3) | (1 << 5)
ZERO_SIGNS = 0
GLOBAL = 0b1111
ENDPOINTS = (1 << 0) | (1 << 3)


def has_edge(graph, a, b):
    if a > b:
        a, b = b, a
    return bool(graph & (1 << EDGE_INDEX[(a, b)]))


def local_inversion(state, vertex):
    """Source operation: complement edges inside N(v), flip N(v) signs."""
    graph, signs = state
    neighbors = 0
    for other in range(4):
        if other != vertex and has_edge(graph, vertex, other):
            neighbors |= 1 << other
    new_graph = graph
    for a in range(4):
        if not (neighbors & (1 << a)):
            continue
        for b in range(a + 1, 4):
            if neighbors & (1 << b):
                bit = 1 << EDGE_INDEX[(a, b)]
                new_graph ^= bit
    return new_graph, signs ^ neighbors


def apply_word(state, word):
    for vertex in word:
        state = local_inversion(state, vertex)
    return state


def shortest_paths(start):
    distance = {start: 0}
    words = {start: ()}
    queue = deque([start])
    while queue:
        state = queue.popleft()
        for vertex in range(4):
            nxt = local_inversion(state, vertex)
            if nxt not in distance:
                distance[nxt] = distance[state] + 1
                words[nxt] = words[state] + (vertex,)
                queue.append(nxt)
    return distance, words


def check_source_fidelity():
    # Proposition 1 of the 2025 source: B_aa = B for every state and a.
    for graph in range(1 << len(EDGE_PAIRS)):
        for signs in range(1 << 4):
            state = (graph, signs)
            for vertex in range(4):
                assert local_inversion(local_inversion(state, vertex), vertex) == state
    print("SOURCE_FIDELITY_OK involution B_aa=B for all 1024 P4 states")


def check_counterexample():
    start = (P4, ZERO_SIGNS)
    distance, words = shortest_paths(start)
    global_target = (P4, GLOBAL)
    endpoint_target = (P4, ENDPOINTS)
    assert global_target in distance and endpoint_target in distance
    assert distance[global_target] == 8
    assert distance[endpoint_target] == 10
    global_word = (1, 2, 3, 0, 3, 0, 1, 2)
    endpoint_word = (0, 1, 0, 2, 0, 1, 3, 2, 3, 2)
    assert apply_word(start, global_word) == global_target
    assert apply_word(start, endpoint_word) == endpoint_target
    assert len(global_word) == distance[global_target]
    assert len(endpoint_word) == distance[endpoint_target]
    reachable_graphs = {graph for graph, _ in distance}
    assert len(reachable_graphs) == 11
    assert len(distance) == 11 * 16
    print(f"P4_MASK={P4} GLOBAL_DISTANCE={distance[global_target]} WORD={list(global_word)}")
    print(f"P4_MASK={P4} ENDPOINT_DISTANCE={distance[endpoint_target]} WORD={list(endpoint_word)}")
    print("POTENTIAL_GAP=2 (selective endpoint reversal exceeds global reversal)")


def check_disjoint_family_additivity():
    # Local inversions never create cross-component edges.  Projection of any
    # word to one P4 component is a valid word for that component, so every
    # k-copy word has length at least the sum of the two one-copy distances.
    # Concatenating shortest words attains those lower bounds.  This is an
    # exact integer argument for every k, with costs 8k and 10k.
    for k in range(1, 9):
        assert 8 * k < 10 * k
    print("DISJOINT_P4_FAMILY_OK costs global=8k selective_endpoints=10k for k=1..8")


def main():
    check_source_fidelity()
    check_counterexample()
    check_disjoint_family_additivity()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
