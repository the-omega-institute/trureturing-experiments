#!/usr/bin/env python3
"""Exact BFS checks for the modified indexed BIBD graph in arXiv:2611.05762.

All arithmetic is integer arithmetic.  The general formula and the independently
implemented BFS are compared at every ordered vertex pair of each tested graph.
This is finite experimental evidence, not a proof for arbitrary designs.
"""

from collections import Counter, deque
from itertools import combinations, product


def parameters(v, blocks):
    assert blocks
    k = len(blocks[0])
    assert v >= k >= 2
    assert all(len(b) == k and len(set(b)) == k for b in blocks)
    assert all(0 <= x < v for b in blocks for x in b)
    replication = [sum(x in b for b in blocks) for x in range(v)]
    pairs = [sum(x in b and y in b for b in blocks)
             for x, y in combinations(range(v), 2)]
    assert len(set(replication)) == len(set(pairs)) == 1
    r, lam = replication[0], pairs[0]
    assert lam >= 1
    assert v * r == len(blocks) * k
    assert lam * (v - 1) == r * (k - 1)
    return len(blocks), r, k, lam


def graph(v, blocks):
    vertices = [("p", x) for x in range(v)]
    vertices += [("f", i, x) for i, b in enumerate(blocks) for x in b]
    adj = {a: set() for a in vertices}

    def edge(a, b):
        assert a != b
        adj[a].add(b)
        adj[b].add(a)

    for i, b in enumerate(blocks):
        for x in b:
            edge(("p", x), ("f", i, x))
        for x, y in combinations(b, 2):
            edge(("f", i, x), ("f", i, y))
    return vertices, adj


def bfs(adj, start):
    distance, count = {start: 0}, {start: 1}
    queue = deque([start])
    while queue:
        a = queue.popleft()
        for b in adj[a]:
            if b not in distance:
                distance[b], count[b] = distance[a] + 1, 0
                queue.append(b)
            if distance[b] == distance[a] + 1:
                count[b] += count[a]
    assert len(distance) == len(adj)
    return {a: (distance[a], count[a]) for a in adj}


def formula(blocks, lam, a, b):
    if a == b:
        return (0, 1), "equal"
    if a[0] == b[0] == "p":
        return (3, lam), "point-point"
    if a[0] == "f" and b[0] == "p":
        return formula(blocks, lam, b, a)
    if a[0] == "p":
        x, i, y = a[1], b[1], b[2]
        if x == y:
            return (1, 1), "point-flag incident label"
        if x in blocks[i]:
            return (2, 1), "point-flag same block"
        return (4, lam), "point-flag outside block"
    i, x, j, y = a[1], a[2], b[1], b[2]
    if i == j:
        return (1, 1), "flags same block"
    if x == y:
        return (2, 1), "flags same label"
    e = int(x in blocks[j]) + int(y in blocks[i])
    s = len(set(blocks[i]) & set(blocks[j]))
    if e:
        return (3, e), "flags cross incidence e=" + str(e)
    if s:
        return (4, s), "flags intersection only"
    return (5, lam), "flags disjoint blocks"


def check(name, v, blocks, verbose=True):
    b, r, k, lam = parameters(v, blocks)
    vertices, adj = graph(v, blocks)
    # Fidelity: source Proposition 6.10.1, equations (8) and (9).
    assert len(vertices) == v + b * k == v * (r + 1)
    assert all(len(adj[("p", x)]) == r for x in range(v))
    assert all(len(adj[("f", i, x)]) == k
               for i, block in enumerate(blocks) for x in block)
    cases, diameter, max_count = Counter(), 0, 0
    for a in vertices:
        observed = bfs(adj, a)
        for z in vertices:
            expected, case = formula(blocks, lam, a, z)
            assert observed[z] == expected, (name, a, z, observed[z], expected)
            cases[case] += 1
            diameter = max(diameter, observed[z][0])
            max_count = max(max_count, observed[z][1])
    if verbose:
        print(f"{name}: (v,b,r,k,lambda)={(v,b,r,k,lam)} "
              f"vertices={len(vertices)} ordered_pairs={len(vertices)**2} "
              f"diameter={diameter} max_sigma={max_count}")
        print("  cases=" + repr(dict(sorted(cases.items()))))
    return len(vertices) ** 2, cases


def small_designs(v_max, multiplicity_max):
    """Every labelled block multiset, modulo permutation of block indices."""
    for v in range(2, v_max + 1):
        for k in range(2, v + 1):
            possible = tuple(combinations(range(v), k))
            for multiplicities in product(range(multiplicity_max + 1), repeat=len(possible)):
                blocks = tuple(b for b, m in zip(possible, multiplicities) for _ in range(m))
                if not blocks:
                    continue
                replication = [sum(x in b for b in blocks) for x in range(v)]
                if len(set(replication)) != 1:
                    continue
                pairs = [sum(x in b and y in b for b in blocks)
                         for x, y in combinations(range(v), 2)]
                if min(pairs) >= 1 and len(set(pairs)) == 1:
                    yield v, blocks


def six_point_design():
    triples = tuple(combinations(range(6), 3))
    pairs = tuple(combinations(range(6), 2))
    for blocks in combinations(triples, 10):
        counts = Counter(pair for b in blocks for pair in combinations(b, 2))
        if all(counts[pair] == 2 for pair in pairs):
            return blocks
    raise AssertionError("No simple 2-(6,3,2) design found")


def main():
    fano = ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5),
            (1, 4, 6), (2, 3, 6), (2, 4, 5))
    six = six_point_design()
    print("Exact integer BFS; the zero-edge path has multiplicity one.")
    print("Source fidelity: Proposition 6.10.1 (vertex numbers and both degrees).")
    print("Constructed 2-(6,3,2) blocks=" + repr(six))
    examples = [("Fano plane 2-(7,3,1)", 7, fano),
                ("simple 2-(6,3,2)", 6, six),
                ("twice indexed Fano blocks", 7, fano + fano),
                ("complete 2-(6,3,4)", 6, tuple(combinations(range(6), 3))),
                ("complete 2-(4,2,1)", 4, tuple(combinations(range(4), 2))),
                ("boundary k=v with three indexed copies", 3, ((0, 1, 2),) * 3)]
    all_cases, total_pairs = Counter(), 0
    for name, v, blocks in examples:
        n, cases = check(name, v, blocks)
        total_pairs += n
        all_cases.update(cases)
    for v_max, multiplicity_max in [(5, 1), (4, 2)]:
        designs = pairs = 0
        for v, blocks in small_designs(v_max, multiplicity_max):
            n, cases = check("exhaustive", v, blocks, verbose=False)
            designs += 1
            pairs += n
            all_cases.update(cases)
        print(f"Exhaustive v=2..{v_max}, each block multiplicity=0..{multiplicity_max}: "
              f"designs={designs} ordered_pairs={pairs}")
        total_pairs += pairs
    required_cases = {"equal", "point-point", "point-flag incident label",
                      "point-flag same block", "point-flag outside block",
                      "flags same block", "flags same label",
                      "flags cross incidence e=1", "flags cross incidence e=2",
                      "flags intersection only", "flags disjoint blocks"}
    assert set(all_cases) == required_cases
    print("All eleven table branches exercised; total ordered-pair checks=" + str(total_pairs))
    print("Finite checks do not establish the universally quantified table or its novelty.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
