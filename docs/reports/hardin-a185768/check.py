#!/usr/bin/env python3
"""Exact source fidelity and seeded transfer-polynomial checks for OEIS A185768.

This is a computational certificate supporting a proposed universal proof.
It is not a Lean proof. Run with Python 3; no external packages or network needed.
"""

import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = json.loads((HERE / 'source.json').read_text())
C = [1, -1, -179, -341, 6316, 6138, -92272, 87872, 479803,
     -1815915, 852208, 6873004, -8175787, -6036792, 5492010,
     12389568, -3504328, -38808492, 4735254, 36745560, 36433267,
     -12744011, -28480781, -9986028, 2174189, 3904760, 845298,
     118044, 60864, 3296, 32]


def literal_valid(rows):
    """Literal source predicate on an actual array of rows."""
    traces = [[rows[r][j] + rows[r + 1][j + 1] for j in range(8)]
              for r in range(len(rows) - 1)]
    return (all(t[j] != t[j + 1] for t in traces for j in range(7))
            and all(traces[r][j] != traces[r + 1][j]
                    for r in range(len(traces) - 1) for j in range(8)))


def tuple_model():
    rows = list(itertools.product((0, 1), repeat=9))
    states = [(a, b) for a in rows for b in rows if literal_valid((a, b))]
    by_first = {a: [] for a in rows}
    for i, (a, _) in enumerate(states):
        by_first[a].append(i)
    edges = [[j for j in by_first[b] if literal_valid((a, b, states[j][1]))]
             for a, b in states]
    return rows, states, edges


def mask_model():
    """Independent trace encoding: disjoint bit sets for trace values 0 and 2."""
    full = (1 << 8) - 1
    def labels(a, b):
        top, bottom = a & full, (b >> 1) & full
        zero = full ^ (top | bottom)
        two = top & bottom
        one = full ^ (zero | two)
        return zero, one, two
    states = []
    lab = []
    for a in range(1 << 9):
        for b in range(1 << 9):
            t = labels(a, b)
            if all((x & (x >> 1)) == 0 for x in t):
                states.append((a, b))
                lab.append(t)
    by_first = [[] for _ in range(1 << 9)]
    for i, (a, _) in enumerate(states):
        by_first[a].append(i)
    incoming = [[] for _ in states]
    edges = [[] for _ in states]
    for i, (_, b) in enumerate(states):
        for j in by_first[b]:
            if all((x & y) == 0 for x, y in zip(lab[i], lab[j])):
                incoming[j].append(i)
                edges[i].append(j)
    return states, edges, incoming


def step(v, edges):
    return [sum(v[j] for j in e) for e in edges]


def residual(k, vectors):
    return [sum(c * vectors[k - j][i] for j, c in enumerate(C))
            for i in range(len(vectors[0]))]


def main():
    assert len(C) == 31 and C[0] == 1
    # Parse the printed formula independently, rather than merely duplicate C.
    import re
    text = SOURCE['formula']
    terms = re.findall(r'([+-]?)(?:(\d+)\*)?a\(n-(\d+)\)',
                       text.split('=', 1)[1])
    parsed = [1] + [None] * 30
    for sign, coefficient, lag in terms:
        parsed[int(lag)] = -(int(coefficient or '1') * (-1 if sign == '-' else 1))
    assert parsed == C
    assert text.endswith('for n>36')
    print('source=A185768; revision=%s; updated=%s' %
          (SOURCE['revision'], SOURCE['updated']))
    print('C=' + str(C))
    print('printed_formula_coefficients_and_threshold=OK')

    rows, states, edges = tuple_model()
    assert len(states) == 11400
    assert sum(map(len, edges)) == 99464
    # Independent literal enumeration of every 18-cell source array (n=1).
    brute = sum(literal_valid((a, b)) for a in rows for b in rows)
    assert brute == SOURCE['displayed_terms'][0] == len(states)
    print('literal_two_row_arrays_checked=%s; valid=%s' % (2 ** 18, brute))

    masks, mask_edges, incoming = mask_model()
    def bits(row):
        return sum(x << j for j, x in enumerate(row))
    where = {s: i for i, s in enumerate(masks)}
    perm = [where[(bits(a), bits(b))] for a, b in states]
    assert len(set(perm)) == len(states) == len(masks)
    for i, es in enumerate(edges):
        assert {perm[j] for j in es} == set(mask_edges[perm[i]])
    print('tuple_and_bitmask_state_edge_bijection=OK')
    print('states=%s; directed_edges=%s' % (len(states), sum(map(len, edges))))
    encoding = json.dumps([[list(a), list(b), es]
                           for (a, b), es in zip(states, edges)],
                          separators=(',', ':')).encode()
    print('canonical_tuple_graph_sha256=' + hashlib.sha256(encoding).hexdigest())

    expected = [int(x) for x in SOURCE['bfile_terms']]
    assert len(expected) == 200
    vectors = []
    totals = []
    backward = [1] * len(states)
    forward = [1] * len(masks)
    for k in range(200):
        if k <= 38:
            vectors.append(backward)
        total = sum(backward)
        assert total == sum(forward) == expected[k], ('source mismatch', k + 1)
        totals.append(total)
        if k < 199:
            backward = step(backward, edges)
            forward = step(forward, incoming)
    assert totals[:14] == SOURCE['displayed_terms']
    print('displayed_prefix_1_to_14=OK')
    print('published_bfile_1_to_200_tuple_backward_and_mask_forward=OK')
    print('first_14=' + str(totals[:14]))

    example_rows = [[int(x) for x in re.findall('[01]', r)]
                    for r in SOURCE['examples']]
    assert all(len(r) == 18 for r in example_rows)
    for start in (0, 9):
        assert literal_valid([r[start:start + 9] for r in example_rows])
    print('two_displayed_3_by_9_examples=OK')

    R = {k: residual(k, vectors) for k in (36, 37, 38)}
    for k in (36, 37, 38):
        print('r_%s: sum=%s; nonzero_coordinates=%s; max_abs=%s' %
              (k, sum(R[k]), sum(x != 0 for x in R[k]), max(map(abs, R[k]))))
    assert sum(R[36]) == sum(R[37]) == 0
    assert all(x == 0 for x in R[38])
    assert any(x != 0 for x in R[36])
    assert step(R[36], edges) == R[37]
    assert step(R[37], edges) == R[38]
    print('seeded_identity_T_to_8_P_of_T_one_zero=OK')
    print('two_scalar_boundary_identities_n_37_and_38=OK')
    print('residual_linearity_r37_T_r36_and_r38_T_r37=OK')
    for n in range(37, 201):
        assert sum(c * totals[n - j - 1] for j, c in enumerate(C)) == 0, n
    print('scalar_recurrence_all_tested_n_37_to_200=OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
