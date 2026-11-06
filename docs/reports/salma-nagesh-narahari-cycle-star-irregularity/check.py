#!/usr/bin/env python3
"""Exact finite checks for Salma–Nagesh–Narahari Conjecture 3.1.

Only integer arithmetic is used. Finite checks do not prove the unbounded claim.
The search solver does not use the proposed path-label formula.
"""


def edges(k, m):
    """Source indices: root 0, leaves 1..m, cycle path m+1..k+m-1."""
    n = k + m
    return ([(0, i) for i in range(1, m + 1)]
            + [(0, m + 1), (0, n - 1)]
            + [(i, i + 1) for i in range(m + 1, n - 1)])


def valid(k, m, bound, labels):
    assert len(labels) == k + m
    assert all(isinstance(x, int) and 1 <= x <= bound for x in labels)
    graph_edges = edges(k, m)
    assert len(graph_edges) == len(set(graph_edges)) == k + m
    weights = [labels[u] + labels[v] for u, v in graph_edges]
    assert len(set(weights)) == k + m, (k, m, bound, labels, weights)
    return weights


def target(k, m):
    return max(m + 2, (k + m + 2) // 2)


def source_value(k, m):
    """The source's Theorems 3.2–3.6 and Conjecture 3.1, literally by cases."""
    n = k + m
    if k == 3:
        return n - 1
    if k == 4:
        return n - 2
    if k == 5:
        return n - 2 if n == 6 else n - 3
    if k == 6:
        return n - 3 if n in (7, 8) else n - 4
    if k == 7:
        return n - 3 if n == 8 else n - 4 if n in (9, 10) else n - 5
    assert k >= 8 and m >= 1
    if k + 1 <= n <= 2 * k - 4:
        return (n + 2) // 2
    assert n >= 2 * k - 3
    return n - k + 2


def construct(k, m):
    assert k >= 4 and m >= 1
    bound = target(k, m)
    if bound >= k - 1:
        path = [(j + 1) // 2 for j in range(1, k)]
        internal = [a + b for a, b in zip(path, path[1:])]
        assert internal == list(range(2, k))
        endpoint = path[-1]
        assert 2 <= endpoint <= bound
        reserve = [x for x in range(1, bound + 1) if x not in (1, endpoint)]
        assert len(reserve) == bound - 2 >= m
        branch = "large"
    else:
        assert bound <= k - 2 and k <= 2 * bound - 2
        path = [(j + 1) // 2 + int(j > bound and (j - bound) % 2 == 1)
                for j in range(1, k)]
        internal = [a + b for a, b in zip(path, path[1:])]
        assert internal == [x for x in range(2, k + 1) if x != bound + 1]
        assert path[0] == 1
        endpoint = path[-1]
        assert k - bound + 1 <= endpoint <= bound
        assert bound + endpoint >= k + 1
        reserve = [x for x in range(k - bound + 1, bound + 1) if x != endpoint]
        assert len(reserve) == 2 * bound - k - 1 >= m
        branch = "small"
    labels = [bound] + reserve[:m] + path
    valid(k, m, bound, labels)
    return branch, labels


def search_labeling(k, m, bound):
    """Complete finite feasibility search, with necessary counting prunes.

    Assign the root and the cycle path by backtracking, pruning repeated sums.
    At the end, leaf labels are exactly the labels whose root sums are unused;
    any m of them suffice. Vertices may share labels. No formula from construct
    is used here. The degree and weight-range prunes are graph invariants.
    """
    graph_edges = edges(k, m)
    degree = [0] * (k + m)
    for u, v in graph_edges:
        degree[u] += 1
        degree[v] += 1
    if bound < max(degree) or len(graph_edges) > 2 * bound - 1:
        return None
    for root in range(bound, 0, -1):
        path = []
        used = set()

        def extend():
            if len(path) == k - 1:
                closing = root + path[-1]
                if closing in used:
                    return None
                occupied = used | {closing}
                leaves = [x for x in range(1, bound + 1) if root + x not in occupied]
                if len(leaves) < m:
                    return None
                return [root] + leaves[:m] + path[:]
            previous = path[-1] if path else root
            for value in range(1, bound + 1):
                weight = previous + value
                if weight in used:
                    continue
                path.append(value)
                used.add(weight)
                found = extend()
                used.remove(weight)
                path.pop()
                if found is not None:
                    return found
            return None

        found = extend()
        if found is not None:
            return found
    return None


def main():
    # Fidelity: reproduce exactly the labeling in the source's Theorem 3.2.
    for m in range(1, 41):
        labels = [1] + list(range(1, m + 3))
        weights = valid(3, m, m + 2, labels)
        assert weights == list(range(2, m + 2)) + [m + 2, m + 3, 2 * m + 3]
    print("FIDELITY Theorem 3.2: k=3, 1<=m<=40, 40 source labelings verified")

    # A second implementation of the source graph: cyclic order 0..k-1,
    # followed by leaves k..k+m-1. Check the explicit relabeling to source indices.
    for k in range(3, 13):
        for m in range(1, 7):
            alternate = [(i, (i + 1) % k) for i in range(k)]
            alternate += [(0, k + i) for i in range(m)]
            mapping = [0] + list(range(m + 1, k + m)) + list(range(1, m + 1))
            mapped = {tuple(sorted((mapping[a], mapping[b]))) for a, b in alternate}
            assert mapped == {tuple(sorted(e)) for e in edges(k, m)}
    print("GRAPH_FIDELITY 60 source/cyclic edge-set isomorphisms verified")

    # Full finite minimum search, independently of the proposed construction.
    count = 0
    for k in range(3, 13):
        for m in range(1, 7):
            expected = source_value(k, m)
            witnesses = [(bound, search_labeling(k, m, bound))
                         for bound in range(1, expected + 1)]
            feasible = [(bound, labels) for bound, labels in witnesses if labels is not None]
            assert feasible and feasible[0][0] == expected
            valid(k, m, *feasible[0])
            count += 1
    print(f"EXACT_MINIMA 3<=k<=12, 1<=m<=6: {count} minima equal source values")

    counts = {"large": 0, "small": 0}
    source_cases = 0
    for k in range(4, 161):
        for m in range(1, 81):
            branch, _ = construct(k, m)
            counts[branch] += 1
            assert target(k, m) == source_value(k, m)
            source_cases += int(k >= 8)
    assert all(counts.values())
    print(f"CONSTRUCTION 4<=k<=160, 1<=m<=80: {sum(counts.values())} labelings verified")
    print(f"BRANCHES large={counts['large']} small={counts['small']}; exact inventories/reserves verified")
    print(f"CONJECTURE_PIECEWISE 8<=k<=160, 1<=m<=80: {source_cases} identities verified")
    print("ARITHMETIC integers only; finite checks are not an unbounded proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
