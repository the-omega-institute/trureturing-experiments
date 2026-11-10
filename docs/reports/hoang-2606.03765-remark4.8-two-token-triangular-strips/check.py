#!/usr/bin/env python3
"""Exact independent checks for Hoang arXiv:2606.03765v1, Remark 4.8.

Only Python's standard library is used. All graph decisions use integers,
sets, exhaustive finite enumeration and exact backtracking. No numerical
approximation or floating point enters a check.
"""
from collections import Counter
from itertools import combinations, permutations


def edge(a, b):
    return (min(a, b), max(a, b))


def graph(order, edges):
    E = set(edges)
    return tuple(tuple(i != j and edge(i, j) in E for j in range(order))
                 for i in range(order))


def triangular_strip(n):
    vertices = [(r, i) for r in range(n) for i in range(3)]
    return tuple(tuple(a != b and (a[0] == b[0] or
                       (a[1] == b[1] and abs(a[0] - b[0]) == 1))
                       for b in vertices) for a in vertices)


def complement(order, edges):
    return set(combinations(range(order), 2)) - set(edges)


def token_sliding(order, host_edges):
    """Literal source definition: ALL independent two-sets, singleton moves."""
    E = set(host_edges)
    states = [frozenset(e) for e in combinations(range(order), 2) if e not in E]
    H = []
    for I in states:
        row = []
        for J in states:
            removed, added = I - J, J - I
            row.append(len(removed) == len(added) == 1 and
                       edge(next(iter(removed)), next(iter(added))) in E)
        H.append(tuple(row))
    return tuple(H)


def root_states(root_edges):
    """Complement-root rule; independent implementation of adjacency."""
    E = set(root_edges)
    states = sorted(E)
    return tuple(tuple(len(set(e) & set(f)) == 1 and
                       tuple(sorted(set(e) ^ set(f))) not in E
                       for f in states) for e in states)


def connected(H):
    if not H:
        return False
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for w, adjacent in enumerate(H[v]):
            if adjacent and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(H)


def isomorphic(A, B):
    if len(A) != len(B):
        return False
    size = len(A)
    da, db = list(map(sum, A)), list(map(sum, B))
    if sorted(da) != sorted(db):
        return False
    order = sorted(range(size), key=lambda i: -da[i])
    assigned, used = {}, set()

    def visit(t):
        if t == size:
            return True
        v = order[t]
        for w in range(size):
            if w in used or da[v] != db[w]:
                continue
            if any(A[v][i] != B[w][j] for i, j in assigned.items()):
                continue
            assigned[v] = w
            used.add(w)
            if visit(t + 1):
                return True
            del assigned[v]
            used.remove(w)
        return False

    return visit(0)


def canonical(H):
    size = len(H)
    pairs = list(combinations(range(size), 2))
    return (size, min(sum(int(H[p[i]][p[j]]) << k
                          for k, (i, j) in enumerate(pairs))
                      for p in permutations(range(size))))


def induced_cycles(H):
    """Enumerate chordless cycles by induced paths, without a layer argument."""
    cycles = set()
    size = len(H)
    for start in range(size):
        def extend(path):
            last = path[-1]
            for w in range(start + 1, size):
                if w in path or not H[last][w]:
                    continue
                if H[start][w]:
                    if len(path) >= 2 and not any(H[u][w] for u in path[1:-1]):
                        c = tuple(path + [w])
                        cycles.add(min(c, (start,) + tuple(reversed(c[1:]))))
                elif not any(H[u][w] for u in path[:-1]):
                    extend(path + [w])
        for second in range(start + 1, size):
            if H[start][second]:
                extend([start, second])
    return Counter(map(len, cycles))


def realization_search(H):
    """Complete search for connected nonempty targets, up to root renaming.

    Target vertices receive distinct root edges. A target adjacency requires
    intersection in one endpoint and forbids the edge on the other endpoints.
    A target nonadjacency with intersecting labels requires that blocking edge.
    At completion, selected root edges must exhaust all required edges.

    Fresh endpoints are introduced in increasing order of first appearance;
    swapping the two fresh endpoints of an edge has no effect. Every connected
    target with m states has connected root support and at most m+1 active root
    vertices. Isolated root vertices, i.e. universal host vertices, never occur
    in a state and can be removed. The search therefore does not impose a host
    order restriction on realizability.

    The sole additional pruning is a necessary bipartite matching condition:
    every required but unselected edge must occupy a distinct unfilled target
    vertex, compatibly with all selected labels. It does not use the proposed
    induced-cycle obstruction or the root-claw/interface classification.
    """
    assert connected(H)
    m = len(H)
    nodes = Counter()

    def visit(labels, required, forbidden, order):
        t = len(labels)
        nodes[t] += 1
        if t == m:
            return labels if required <= set(labels) else None
        pending = required - set(labels)
        slots = []
        for e in pending:
            possible = []
            for q in range(t, m):
                compatible = True
                for j, f in enumerate(labels):
                    common = set(e) & set(f)
                    if H[q][j]:
                        if len(common) != 1:
                            compatible = False
                            break
                        if tuple(sorted(set(e) ^ set(f))) in required:
                            compatible = False
                            break
                    elif len(common) == 1:
                        if tuple(sorted(set(e) ^ set(f))) in forbidden:
                            compatible = False
                            break
                if compatible:
                    possible.append(q)
            if not possible:
                return None
            slots.append(possible)
        matching = {}

        def augment(k, seen):
            for q in slots[k]:
                if q in seen:
                    continue
                seen.add(q)
                if q not in matching or augment(matching[q], seen):
                    matching[q] = k
                    return True
            return False

        if any(not augment(k, set())
               for k in sorted(range(len(slots)), key=lambda k: len(slots[k]))):
            return None
        for e in combinations(range(order + 2), 2):
            fresh = sorted(set(e) - set(range(order)))
            if fresh and fresh != list(range(order, order + len(fresh))):
                continue
            if e in forbidden or e in labels:
                continue
            nr, nf = required | {e}, set(forbidden)
            compatible = True
            for j, f in enumerate(labels):
                common = set(e) & set(f)
                if H[t][j]:
                    if len(common) != 1:
                        compatible = False
                        break
                    nf.add(tuple(sorted(set(e) ^ set(f))))
                elif len(common) == 1:
                    nr.add(tuple(sorted(set(e) ^ set(f))))
            if not compatible or nr & nf or len(nr) > m:
                continue
            answer = visit(labels + (e,), nr, nf, max(order, max(e) + 1))
            if answer is not None:
                return answer
        return None

    result = visit(((0, 1),), {(0, 1)}, set(), 2)
    if result is not None:
        # Recompute every adjacency from the final, exhaustive root edge set.
        E = set(result)
        actual = tuple(tuple(len(set(e) & set(f)) == 1 and
                             tuple(sorted(set(e) ^ set(f))) not in E
                             for f in result) for e in result)
        assert actual == H
    return result, sum(nodes.values())


def main():
    # Source Corollary 4.4: K3 disjoint K1 and K3 disjoint K2.
    clique3 = {(0, 1), (0, 2), (1, 2)}
    for n, order, E in [(1, 4, clique3), (2, 5, clique3 | {(3, 4)})]:
        assert isomorphic(token_sliding(order, E), triangular_strip(n))
        print(f'FIDELITY Corollary 4.4: K3 disjoint K{n} realizes H_{n}')

    # Independent brute-force enumeration of every labelled host up to order 6.
    small_realizable = set()
    hits = Counter()
    host_count = 0
    for order in range(7):
        pairs = list(combinations(range(order), 2))
        for mask in range(1 << len(pairs)):
            root = {e for k, e in enumerate(pairs) if (mask >> k) & 1}
            H = root_states(root)
            host_count += 1
            if order <= 5:
                assert H == token_sliding(order, complement(order, root))
                if 1 <= len(H) <= 4 and connected(H):
                    small_realizable.add(canonical(H))
            if len(H) in (3, 6, 9, 12, 15):
                n = len(H) // 3
                if isomorphic(H, triangular_strip(n)):
                    hits[n] += 1
    assert hits[1] > 0 and hits[2] > 0
    assert all(hits[n] == 0 for n in (3, 4, 5))
    print(f'BRUTE_FORCE all labelled hosts of orders 0..6: {host_count}; '
          f'H_n realization counts {dict(sorted(hits.items()))}')
    print('FIDELITY complement/state adjacency agrees for every host of orders 0..5')

    # Every connected target of at most 4 vertices has active root order <= 5.
    targets = {}
    for size in range(1, 5):
        pairs = list(combinations(range(size), 2))
        for mask in range(1 << len(pairs)):
            H = graph(size, {e for k, e in enumerate(pairs) if (mask >> k) & 1})
            if connected(H):
                targets[canonical(H)] = H
    for key, H in sorted(targets.items()):
        answer, _ = realization_search(H)
        assert (answer is not None) == (key in small_realizable)
    print(f'SEARCH_CROSS_CHECK {len(targets)} connected target isomorphism classes '
          'of orders 1..4 agree with exhaustive hosts of orders 0..5')

    # Independent direct chordless-cycle enumeration, with a positive control.
    assert induced_cycles(graph(6, {edge(i, (i + 1) % 6) for i in range(6)})) == {6: 1}
    for n in range(1, 11):
        counts = induced_cycles(triangular_strip(n))
        expected = Counter({3: n})
        if n >= 2:
            expected[4] = 3 * (n - 1)
        assert counts == expected
        print(f'INDUCED_CYCLES H_{n}: {dict(sorted(counts.items()))}; none of length >=5')

    # Six blockers forced by a same-center interface. Other root edges outside
    # these seven vertices cannot change adjacency among these six states.
    xs, ys = (1, 2, 3), (4, 5, 6)
    interface_root = {edge(0, v) for v in xs + ys}
    interface_root |= {edge(xs[i], ys[j]) for i in range(3) for j in range(3) if i != j}
    cycle = [(1, 5), (3, 5), (3, 4), (2, 4), (2, 6), (1, 6)]
    states = sorted(interface_root)
    H = root_states(interface_root)
    indices = [states.index(e) for e in cycle]
    induced = tuple(tuple(H[i][j] for j in indices) for i in indices)
    assert induced_cycles(induced) == {6: 1}
    print('LOCAL_OBSTRUCTION six forced blockers induce exactly C6')

    # Complete small-n search, not a Lean proof of the infinite statement.
    for n in range(1, 8):
        answer, nodes = realization_search(triangular_strip(n))
        assert (answer is not None) == (n <= 2)
        print(f'EXACT_REALIZATION_SEARCH n={n}; states={3*n}; '
              f'active_root_order_bound={3*n+1}; search_nodes={nodes}; '
              f'realizable={answer is not None}')
    print('LIMIT: finite checks do not prove the uniform claim for every n>=6')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
