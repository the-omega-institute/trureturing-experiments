#!/usr/bin/env python3
"""Exact checks for Turau's first radius conjecture, arXiv:2002.10752v4.

Standard library only. Graphs are symmetric lists of sets of integer vertices.
An independent bit-mask implementation checks the literal directed-edge update.
This is experimental evidence, not a Lean proof or an exhaustive literature claim.
"""
from collections import deque
from itertools import combinations


def graph(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        adj[u].add(v)
        adj[v].add(u)
    return adj


def distances(adj, sources):
    d = [None] * len(adj)
    q = deque(sources)
    for v in sources:
        d[v] = 0
    while q:
        u = q.popleft()
        for v in adj[u]:
            if d[v] is None:
                d[v] = d[u] + 1
                q.append(v)
    assert all(x is not None for x in d), 'disconnected graph'
    return d


def radius(adj):
    return min(max(distances(adj, [v])) for v in range(len(adj)))


def is_bipartite(adj):
    color = {}
    for s in range(len(adj)):
        if s in color:
            continue
        color[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    q.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def initial(adj, sources):
    return frozenset((u, v) for u in sources for v in adj[u])


def step(adj, messages):
    incoming = [set() for _ in adj]
    for u, v in messages:
        incoming[v].add(u)
    return frozenset((v, w) for v, heard in enumerate(incoming)
                     if heard for w in adj[v] - heard)


def execution(adj, sources):
    m = initial(adj, sources)
    rounds, seen = [], set()
    while m:
        assert m not in seen, 'nonterminating periodic execution'
        seen.add(m)
        rounds.append(m)
        m = step(adj, m)
    return rounds


def mask_step(neighbors, incoming):
    outgoing = [0] * len(neighbors)
    for v, heard in enumerate(incoming):
        if heard:
            targets = neighbors[v] & ~heard
            while targets:
                bit = targets & -targets
                targets -= bit
                outgoing[bit.bit_length() - 1] |= 1 << v
    return tuple(outgoing)


def mask_initial(adj, sources):
    incoming = [0] * len(adj)
    for s in sources:
        for v in adj[s]:
            incoming[v] |= 1 << s
    return tuple(incoming)


def decode_masks(masks):
    return frozenset((u, v) for v, mask in enumerate(masks)
                     for u in range(len(masks)) if mask & (1 << u))


def fidelity():
    # Source: immediately after Definition 1, complete-graph values.
    instances = 0
    for n in range(3, 10):
        adj = graph(n, combinations(range(n), 2))
        for k in range(1, n + 1):
            expected = 3 if k == 1 else (1 if k == n else 2)
            times = {len(execution(adj, s)) for s in combinations(range(n), k)}
            assert times == {expected}, (n, k, times)
            instances += sum(1 for _ in combinations(range(n), k))
    print(f'FIDELITY complete graphs K_3..K_9: {instances} nonempty source sets; Flood_1=3, Flood_k=2 (1<k<n), Flood_n=1')
    # Source: the bipartite single-source corollary and the odd-cycle example.
    for n in range(2, 17):
        adj = graph(n, ((i, i + 1) for i in range(n - 1)))
        for s in range(n):
            assert len(execution(adj, [s])) == max(s, n - 1 - s)
    for n in range(3, 18, 2):
        adj = graph(n, ((i, (i + 1) % n) for i in range(n)))
        assert len(execution(adj, [0])) == n
    print('FIDELITY paths P_2..P_16: single-source time=eccentricity; odd cycles C_3..C_17: single-source time=n')


def decorated_path(k):
    assert k >= 1
    edges = [(i, i + 1) for i in range(5 * k - 1)]
    for j in range(k):
        l, c, r = 5 * j + 1, 5 * j + 2, 5 * j + 3
        a, b = 5 * k + 2 * j, 5 * k + 2 * j + 1
        edges.extend([(a, l), (a, c), (b, c), (b, r)])
    return graph(7 * k, edges), {5 * j + 2 for j in range(k)}


def predicted_rounds(k):
    m1, m2, m3 = set(), set(), set()
    for j in range(k):
        L, l, c, r, R = (5 * j + t for t in range(5))
        a, b = 5 * k + 2 * j, 5 * k + 2 * j + 1
        m1.update((c, x) for x in (l, a, r, b))
        m2.update([(l, L), (l, a), (a, l), (r, R), (r, b), (b, r)])
        m3.update([(l, L), (l, c), (a, c), (r, R), (r, c), (b, c)])
        if j + 1 < k:
            m3.update([(R, L + 5), (L + 5, R)])
    return [frozenset(m1), frozenset(m2), frozenset(m3)]


def family_checks():
    for k in range(1, 33):
        adj, sources = decorated_path(k)
        n, edges = len(adj), sum(map(len, adj)) // 2
        assert n == 7 * k and edges == 9 * k - 1
        assert max(map(len, adj)) == 4 and not is_bipartite(adj)
        all_d = [distances(adj, [v]) for v in range(n)]
        rad = min(map(max, all_d))
        assert rad == (5 * k) // 2
        # Full path-distance preservation; includes both terminal distances.
        for u in range(5 * k):
            for v in range(5 * k):
                assert all_d[u][v] == abs(u - v)
        assert max(distances(adj, sources)) == 2
        assert k * (1 + max(map(len, adj))) < n  # no k radius-one balls cover
        rounds = execution(adj, sources)
        assert rounds == predicted_rounds(k)
        assert list(map(len, rounds)) == [4 * k, 6 * k, 8 * k - 2]
        assert sum(map(len, rounds)) == 2 * edges
        neighbors = [sum(1 << v for v in ns) for ns in adj]
        masks = mask_initial(adj, sources)
        for m in rounds:
            assert decode_masks(masks) == m
            masks = mask_step(neighbors, masks)
        assert not any(masks)
        gap = rad + k - 1 - 3 * k
        assert gap == k // 2 - 1
        if k >= 4:
            assert gap > 0  # tau=3 gives Flood_k <= 3, enough for refutation
        if k % 2 == 0:
            assert 2 * (3 * k) == 2 * rad + k
        # Source Theorem 4(3), multiplied by 2k, forces integral Flood_k >= 3.
        assert 4 * k < 2 * rad + k <= 6 * k
        if k in (1, 2, 3, 4, 5, 8, 16, 32):
            print(f'FAMILY k={k} vertices={n} edges={edges} radius={rad} r_k=2 witnessed_tau=3 counts={[len(m) for m in rounds]} conjecture_gap={gap}')
    print('FAMILY exact BFS and two independent transmission updates: k=1..32; all tests passed')


def exhaustive_optimality():
    total = 0
    for k in range(1, 5):
        adj, sources = decorated_path(k)
        neighbors = [sum(1 << v for v in ns) for ns in adj]
        count = 0
        # M_3 nonempty for every k-set rules out termination in <= 2 rounds.
        for s in combinations(range(len(adj)), k):
            m1 = mask_initial(adj, s)
            m3 = mask_step(neighbors, mask_step(neighbors, m1))
            assert any(m3), ('two-round source set found', k, s)
            count += 1
        assert len(execution(adj, sources)) == 3
        total += count
        print(f'OPTIMALITY k={k}: all {count} source sets have M_3 nonempty; Flood_k=3')
    print(f'OPTIMALITY total source sets checked={total}; no appeal to the published lower bound in these finite checks')


def main():
    fidelity()
    family_checks()
    exhaustive_optimality()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
