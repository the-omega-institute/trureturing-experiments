#!/usr/bin/env python3
"""Exact checks for the broad v3 ornament definition, arXiv:2606.16001.

All arithmetic is Python integer arithmetic. Exhaustion is a finite check of
an unbounded construction, not a proof of the published existence question.
Run: python3 check.py
"""
from itertools import combinations, permutations, product
import json


def weights(order, edges, labels):
    ans = [0] * order
    for tail, head in edges:
        assert 0 <= tail < order and 0 <= head < order and tail != head
        ans[head] += labels[tail]
        ans[tail] -= labels[head]
    return ans


WHEEL = [(4, 0), (1, 4), (2, 4), (4, 3),
         (0, 1), (0, 2), (3, 1), (3, 2)]


def fidelity():
    # Figure 1(a), transcribed from its TikZ edges, labels 1,...,5.
    assert weights(5, WHEEL, list(range(1, 6))) == [0] * 5
    count = 0
    # Remark 1: a directed cycle has zero total weight and every individual
    # weight is nonzero for distinct labels. Both offsets preserve this.
    for order in range(3, 8):
        edges = [(i, (i + 1) % order) for i in range(order)]
        for p in permutations(range(1, order + 1)):
            for offset in (0, 13):
                w = weights(order, edges, [offset + x for x in p])
                assert sum(w) == 0 and all(x != 0 for x in w)
                count += 1
    return {"figure_1a_wheel_weights": [0] * 5,
            "remark_1_cycle_labelings": count}


def generators(size):
    result = []
    for t in range(size - 4):
        result.extend((
            {t: 1, t + 1: -1, t + 2: -1, t + 3: 1},
            {t: 1, t + 1: -1, t + 3: -1, t + 4: 1},
            {t + 1: 1, t + 2: -1, t + 3: -1, t + 4: 1},
        ))
    for row in result:
        assert set(row.values()) <= {-1, 1}
        assert sum(row.values()) == 0
        assert sum((i + 1) * x for i, x in row.items()) == 0
    assert set().union(*(set(r) for r in result)) == set(range(size))
    return result


def decompose(z, rows):
    size = len(z)
    assert sum(z) == sum((i + 1) * x for i, x in enumerate(z)) == 0
    coeffs = [0] * len(rows)
    for i in range(size - 2):
        q = sum((i + 1 - j) * z[j] for j in range(i + 1))
        t = min(i, size - 5)
        # D_t=R2-R3; D_(t+1)=R1+R3-R2; D_(t+2)=R2-R1.
        multipliers = ((0, 1, -1), (1, -1, 1), (-1, 1, 0))[i - t]
        for j, c in enumerate(multipliers):
            coeffs[3 * t + j] += q * c
    rebuilt = [0] * size
    for c, row in zip(coeffs, rows):
        for i, value in row.items():
            rebuilt[i] += c * value
    assert rebuilt == z
    return coeffs


def six_copies(order, edges):
    size = 6 * order
    # Relabel the six oriented copies by complementary pairs, then use the
    # labels themselves as the vertex coordinates (zero based internally).
    labelings = []
    for q in range(3):
        labelings.extend(([2 * q * order + i + 1 for i in range(order)],
                          [2 * q * order + 2 * order - i for i in range(order)]))
    assert sorted(x for p in labelings for x in p) == list(range(1, size + 1))
    copied_edges = [(p[u] - 1, p[v] - 1) for p in labelings for u, v in edges]
    assert len(set(copied_edges)) == 6 * len(edges)
    assert all((v, u) not in set(copied_edges) for u, v in copied_edges)
    d = weights(size, copied_edges, [1] * size)
    a = weights(size, copied_edges, list(range(1, size + 1)))
    assert sum(d) == sum(a) == 0
    assert sum((i + 1) * d[i] for i in range(size)) == 0
    assert sum((i + 1) * a[i] for i in range(size)) == 0
    rows = generators(size)
    return copied_edges, rows, decompose(d, rows), decompose(a, rows)


def construct(order, edges):
    copied_edges, templates, alpha, beta = six_copies(order, edges)
    size = 6 * order
    offset_count = sum(abs(x) for x in alpha)
    low_rows = {}
    high_rows = []
    delta = beta.copy()
    offset = 0
    for j, coefficient in enumerate(alpha):
        sign = (coefficient > 0) - (coefficient < 0)
        for _ in range(abs(coefficient)):
            high_rows.append((offset, j, sign))
            delta[j] += sign * offset
            offset += 1
    assert offset == offset_count
    maximum = 0
    for j, value in enumerate(delta):
        if value:
            x = maximum + 1 + max(value, 0)
            y = maximum + 1 + max(-value, 0)
            assert x != y and x not in low_rows and y not in low_rows
            low_rows[x], low_rows[y] = (j, 1), (j, -1)
            maximum = max(x, y)
    # A zero-contribution triple for each generator guarantees attachment
    # even for an edgeless H or when all affine coefficients vanish.
    for j in range(len(templates)):
        p, q = maximum + 1, maximum + 2
        r = p + q
        assert len({p, q, r}) == 3
        assert all(x not in low_rows for x in (p, q, r))
        low_rows[p], low_rows[q], low_rows[r] = (j, 1), (j, 1), (j, -1)
        maximum = r
    threshold = max(5, maximum + offset_count + 1)
    cross_edge_count = 0
    for n in (threshold, threshold + 1, threshold + 17):
        selected = low_rows.copy()
        for t, j, sign in high_rows:
            label = n - t
            assert label not in selected
            selected[label] = (j, sign)
        assert all(1 <= label <= n for label in selected)
        # Direct actual edge-weight computation, with the nontrivial old
        # wheel 1,...,5 and isolated old vertices at all other labels.
        old_weights = dict(enumerate(weights(5, WHEEL, list(range(1, 6))), 1))
        new_weights = weights(size, copied_edges, list(range(n + 1, n + size + 1)))
        covered = set()
        transpose_product = [0] * size
        for label, (j, sign) in selected.items():
            row = {i: sign * c for i, c in templates[j].items()}
            assert sum(row.values()) == 0
            assert sum((i + 1) * c for i, c in row.items()) == 0
            assert set(row.values()) <= {-1, 1}
            for i, c in row.items():
                # c=1: new -> old, c=-1: old -> new.
                old_weights[label] = old_weights.get(label, 0) + c * (n + i + 1)
                new_weights[i] -= c * label
                transpose_product[i] += c * label
                covered.add(i)
        assert covered == set(range(size))
        assert all(w == 0 for w in old_weights.values())
        assert new_weights == [0] * size
        assert transpose_product == weights(size, copied_edges,
                                            list(range(n + 1, n + size + 1)))
        # All unselected old rows are exactly zero, without iterating to n.
        cross_edge_count = 4 * len(selected)
    return {"old_order_threshold": threshold, "copies": 6,
            "added_vertices": size, "selected_old_labels": len(selected),
            "cross_edges": cross_edge_count, "tested_orders": 3}


def main():
    source = fidelity()
    counts = {}
    maximum_threshold = 0
    for order in range(1, 5):
        pairs = list(combinations(range(order), 2))
        count = 0
        for directions in product((0, 1, -1), repeat=len(pairs)):
            edges = [(u, v) if d == 1 else (v, u)
                     for (u, v), d in zip(pairs, directions) if d]
            result = construct(order, edges)
            maximum_threshold = max(maximum_threshold, result['old_order_threshold'])
            count += 1
        counts[order] = count
    named = {}
    for name, order, edges in (
        ('isolated_vertex', 1, []),
        ('directed_edge', 2, [(0, 1)]),
        ('directed_triangle', 3, [(0, 1), (1, 2), (2, 0)]),
        ('directed_four_cycle', 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
        ('inward_K_1_4', 5, [(i, 0) for i in range(1, 5)]),
        ('transitive_tournament_5', 5, list(combinations(range(5), 2))),
    ):
        named[name] = construct(order, edges)
    print(json.dumps({"source_fidelity": source, "all_oriented_graph_counts": counts,
                      "exhaustive_graphs": sum(counts.values()),
                      "maximum_exhaustive_threshold": maximum_threshold,
                      "named_cases": named,
                      "arithmetic": "Python integers only; no floating point",
                      "scope": "Every labeled oriented graph on 1 through 4 vertices; three old orders per graph. Sparse zero rows cover all unused old vertices. These finite checks do not prove the unbounded theorem."},
                     indent=2, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
