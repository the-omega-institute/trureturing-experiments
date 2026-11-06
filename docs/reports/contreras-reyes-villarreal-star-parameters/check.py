#!/usr/bin/env python3
"""Exact source fidelity and star-parameter obstruction checks; Python stdlib only."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
import json


def graph(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        adj[u].add(v)
        adj[v].add(u)
    return adj


def edges(adj):
    return [(u, v) for u in range(len(adj)) for v in sorted(adj[u]) if u < v]


def complement(adj):
    return graph(len(adj), ((u, v) for u, v in combinations(range(len(adj)), 2)
                            if v not in adj[u]))


def connected(adj):
    seen, stack = {0}, [0]
    while stack:
        u = stack.pop()
        for v in adj[u] - seen:
            seen.add(v)
            stack.append(v)
    return len(seen) == len(adj)


def triangles(adj):
    return sum(len(adj[u] & adj[v] & set(range(v + 1, len(adj))))
               for u, v in edges(adj))


def colourings(adj, k, fixed=None, limit=None):
    """Exhaustive DSATUR with no symmetry pruning; optional labelled fixed colours."""
    n = len(adj)
    col = [-1] * n
    for v, c in (fixed or {}).items():
        assert 0 <= c < k
        col[v] = c
    assert all(col[u] < 0 or col[v] < 0 or col[u] != col[v] for u, v in edges(adj))
    answer = []

    def visit():
        if limit is not None and len(answer) >= limit:
            return
        left = [v for v in range(n) if col[v] < 0]
        if not left:
            answer.append(tuple(col))
            return
        v = max(left, key=lambda w: (len({col[u] for u in adj[w] if col[u] >= 0}),
                                    len(adj[w]), -w))
        banned = {col[u] for u in adj[v] if col[u] >= 0}
        for c in range(k):
            if c not in banned:
                col[v] = c
                visit()
        col[v] = -1

    visit()
    return answer


def mycielski_c5():
    es = [(i, (i + 1) % 5) for i in range(5)]
    es += [(5 + i, (i + delta) % 5) for i in range(5) for delta in [-1, 1]]
    es += [(10, 5 + i) for i in range(5)]
    return graph(11, es)


def categorical_k3(adj):
    n = len(adj)
    return graph(3 * n, ((i * n + b, j * n + c)
                        for i in range(3) for j in range(3) if i != j
                        for b, c in edges(adj)
                        for b, c in [(b, c), (c, b)]))


def blowup(adj, t):
    return graph(len(adj) * t, ((u * t + a, v * t + b)
                              for u, v in edges(adj) for a in range(t) for b in range(t)))


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def no_possible_star_rows(adj_x, cols):
    """A star row needs centre coefficient 1 and zero on every complement-neighbour."""
    return all(not (cols[v][r] == 1 and all(cols[w][r] == 0 for w in adj_x[v]))
               for v in range(len(adj_x)) for r in range(2))


def main():
    b = mycielski_c5()
    assert len(edges(b)) == 20 and triangles(b) == 0 and connected(b)
    assert not colourings(b, 3) and colourings(b, 4, limit=1)

    # Source Example 5.10: the 35-edge complement, including its displayed rational forms.
    gb = complement(b)
    source_edges_one_based = [
        (1, 3), (1, 4), (2, 4), (2, 5), (3, 5), (1, 6), (1, 8), (1, 9),
        (2, 7), (2, 9), (2, 10), (3, 6), (3, 8), (3, 10), (4, 6), (4, 7), (4, 9),
        (5, 7), (5, 8), (5, 10), (1, 11), (2, 11), (3, 11), (4, 11), (5, 11),
        (6, 7), (6, 8), (6, 9), (6, 10), (7, 8), (7, 9), (7, 10), (8, 9), (8, 10), (9, 10)]
    assert edges(gb) == sorted((u - 1, v - 1) for u, v in source_edges_one_based)
    independent = [s for size in range(12) for s in combinations(range(11), size)
                   if all(v not in gb[u] for u, v in combinations(s, 2))]
    assert max(map(len, independent)) == 2
    assert not colourings(gb, 5) and colourings(gb, 6, limit=1)
    f1 = list(map(Q, ['5/3', '8/3', '2/9', '10/3', '7/10', '1/3', '10/9', '6/5', '9/8', '0', '1']))
    f2 = list(map(Q, ['4/5', '1/8', '4/7', '1', '4/7', '6/5', '3/10', '4/3', '3/5', '1', '0']))
    rational_cols = list(zip(f1, f2))
    rational_minors = [det(rational_cols[u], rational_cols[v]) for u, v in edges(b)]
    assert all(rational_minors) and all(any(c) for c in rational_cols)

    # All fibre types in the product-rigidity argument, exhaustively and exactly.
    functions = list(product(range(3), repeat=3))
    compatible = [(h, g) for h in functions for g in functions
                  if all(h[i] != g[j] for i in range(3) for j in range(3) if i != j)]
    assert len(functions) == 27 and len(compatible) == 66
    for h, g in compatible:
        if len(set(h)) == 3:
            assert g == h
        else:
            assert len(set(g)) < 3
            assert max(range(3), key=h.count) != max(range(3), key=g.count)

    x = categorical_k3(b)
    g = complement(x)
    assert len(x) == 33 and len(edges(x)) == 120 and triangles(x) == 0
    assert connected(x) and connected(g) and len(edges(g)) == 408
    assert min(map(len, g)) == 22 and max(map(len, g)) == 26
    u, v = edges(x)[0]
    normalized = colourings(x, 3, fixed={u: 0, v: 1})
    assert len(normalized) == 1
    all_cols = {tuple(p[c] for c in normalized[0]) for p in permutations(range(3))}
    expected_cols = {tuple(p[w // 11] for w in range(33)) for p in permutations(range(3))}
    assert all_cols == expected_cols and len(all_cols) == 6
    binary = [(1, 0), (0, 1), (1, 1)]
    parameter_cols = [binary[w // 11] for w in range(33)]
    minors = [det(parameter_cols[a], parameter_cols[c]) for a, c in edges(x)]
    assert set(minors) <= {-1, 1}
    assert all(len({colour[w] for w in x[z]}) == 2
               for colour in all_cols for z in range(33))
    for colour in all_cols:
        assert no_possible_star_rows(x, [binary[c] for c in colour])

    # Extra twins are forced by representative neighbours in both other colour classes.
    family = []
    for t in [1, 2, 3, 4]:
        xt = blowup(x, t)
        gt = complement(xt)
        cols = [parameter_cols[v // t] for v in range(33 * t)]
        assert triangles(xt) == 0 and len(edges(xt)) == 120 * t * t
        assert len(edges(gt)) == (849 * t * t - 33 * t) // 2 and connected(gt)
        assert all(det(cols[a], cols[c]) in [-1, 1] for a, c in edges(xt))
        for colour in all_cols:
            for z in range(len(xt)):
                representative_colours = {colour[w // t] for w in xt[z] if w % t == 0}
                assert representative_colours == set(range(3)) - {colour[z // t]}
            assert no_possible_star_rows(xt, [binary[colour[z // t]] for z in range(len(xt))])
        family.append({'t': t, 'vertices': len(xt), 'X_edges': len(edges(xt)),
                       'G_edges': len(edges(gt)), 'independence_number_G': 2,
                       'independent_ordered_centre_pairs_excluded': 2 * len(edges(xt))})

    report = {
        'source_fidelity': {'source': 'arXiv:2609.30081v2, Example 5.10',
                            'complement_Myc_vertices': 11, 'edges': 35, 'independence_number': 2,
                            'chromatic_number': 6, 'displayed_rational_pair_minors_nonzero': len(rational_minors),
                            'rational_minors': [str(q) for q in rational_minors]},
        'base': {'vertices': 11, 'edges': 20, 'triangles': 0, 'chromatic_number': 4},
        'rigidity_local': {'fibre_functions': 27, 'compatible_ordered_pairs': 66},
        'counterexample': {'vertices': 33, 'X_edges': 120, 'G_edges': 408,
                           'G_min_degree': 22, 'G_max_degree': 26, 'independence_number_G': 2,
                           'normalized_three_colourings': len(normalized),
                           'all_labelled_three_colourings': len(all_cols),
                           'binary_pair_minors': {str(a): minors.count(a) for a in sorted(set(minors))},
                           'possible_star_rows_across_all_colourings': 0},
        'family_checks': family,
        'scope': 'Exact finite corroboration; infinite-family and arbitrary-field arguments are in the preregistration, not kernel-verified here.'}
    print(json.dumps(report, indent=2, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
