#!/usr/bin/env python3
"""Exact finite checks for Ciobanu arXiv:2609.15879, Question 12.3.

These checks validate the proposed finite invariant and source fidelity.
They do not certify the classification of all nominal atoms or categorical
non-equivalence. Only Python integers and finite combinatorial objects occur.
"""
from collections import Counter
from itertools import product


def rgs(n):
    """All set partitions of n labelled slots, as restricted growth strings."""
    def extend(xs, high):
        if len(xs) == n:
            yield tuple(xs)
            return
        for a in range(high + 2):
            yield from extend(xs + [a], max(high, a))
    yield from extend([], -1)


def normalize(xs):
    names = {}
    return tuple(names.setdefault(x, len(names)) for x in xs)


def edge_key(xs):
    """Quotient by vertex renaming and reversal inside each unordered edge."""
    n = len(xs) // 2
    return min(normalize(tuple(xs[2*i + (side ^ flip[i])]
                               for i in range(n) for side in range(2)))
               for flip in product(range(2), repeat=n))


def nominal_keys(n):
    return {edge_key(xs) for xs in rgs(2*n)
            if all(xs[2*i] != xs[2*i+1] for i in range(n))}


def venn_keys(n):
    """Independent classification by multiplicities of vertex membership types."""
    masks = range(1, 1 << n)
    answer = set()
    for counts in product(range(3), repeat=(1 << n)-1):
        if all(sum(c for m, c in zip(masks, counts) if (m >> i) & 1) == 2
               for i in range(n)):
            edges = [[] for _ in range(n)]
            vertex = 0
            for m, c in zip(masks, counts):
                for _ in range(c):
                    for i in range(n):
                        if (m >> i) & 1:
                            edges[i].append(vertex)
                    vertex += 1
            answer.add(edge_key(tuple(a for e in edges for a in e)))
    return answer


def graph_tuple_count(n):
    """Equality partition plus one adjacency bit for every distinct pair."""
    by_distinct_vertices = Counter()
    for xs in rgs(n):
        k = len(set(xs))
        by_distinct_vertices[k] += 1 << (k*(k-1)//2)
    return sum(by_distinct_vertices.values()), dict(sorted(by_distinct_vertices.items()))


def rado_adj(a, b):
    lo, hi = sorted((a, b))
    return a != b and ((hi >> lo) & 1) == 1


def check_source_fidelity():
    # Source Proposition "Sort counts", item (iii): 2^k fresh adjacency types.
    # The high bit forces the witness outside the finite context. Its lower
    # bits encode the requested adjacency vector literally in the Rado model.
    for k in range(9):
        for context in [tuple(range(k)), tuple(2*i+1 for i in range(k))]:
            high = max(context, default=-1) + 1
            patterns = set()
            for mask in range(1 << k):
                witness = (1 << high) + sum(1 << a for i, a in enumerate(context)
                                           if (mask >> i) & 1)
                assert witness not in context
                observed = tuple(int(rado_adj(a, witness)) for a in context)
                assert observed == tuple((mask >> i) & 1 for i in range(k))
                patterns.add(observed)
            assert len(patterns) == 1 << k
        print(f"FIDELITY random-graph fresh sorts k={k}: {1 << k}")


def main():
    check_source_fidelity()
    nominal_expected = [1, 1, 3, 16]
    for n, expected in enumerate(nominal_expected):
        a, b = nominal_keys(n), venn_keys(n)
        assert a == b and len(a) == expected
        print(f"NOM unordered-two-sets power={n}: {len(a)} (partition and Venn checks agree)")
    graph_expected = [1, 1, 3, 15, 127, 1895]
    for n, expected in enumerate(graph_expected):
        total, profile = graph_tuple_count(n)
        assert total == expected
        print(f"RADO vertex power={n}: {total}; by distinct vertices={profile}")
    # Independent finite check of the support-size lower-bound witnesses.
    for k in range(6):
        left = set(range(k))
        right = [set(range(j)) | set(range(k, 2*k-j)) for j in range(k+1)]
        intersections = [len(left & s) for s in right]
        assert intersections == list(range(k+1))
        assert all(len(s) == k for s in right)
        print(f"SUPPORT size={k}: {k+1} distinct square-orbit witnesses")
    ordered_square = {xs for xs in rgs(4) if xs[0] != xs[1] and xs[2] != xs[3]}
    assert len(ordered_square) == 7
    print(f"NOM ordered-distinct-two-tuples square: {len(ordered_square)}")
    assert len(nominal_keys(2)) == graph_tuple_count(2)[0] == 3
    assert len(nominal_keys(3)) == 16 != graph_tuple_count(3)[0] == 15
    print("SEPARATION proposed finite signatures: nominal=(3,16), Rado=(3,15)")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
