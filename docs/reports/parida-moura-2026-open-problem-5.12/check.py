#!/usr/bin/env python3
"""Exact corroboration of Parida--Moura OP 5.12, not a universal proof.

Uses only Python's standard library. Columns are labelled, so an ECFF may
repeat a column within its matching edge. Source fidelity is Proposition
5.11 of arXiv:2605.12634v3. All assertions use integer or finite-set operations.
"""
from itertools import combinations
from math import comb
import json


def middle(n):
    return comb(n, n // 2)


def mask(points):
    return sum(1 << x for x in points)


def check_masks(edges, rows, sperner):
    columns = [column for pair in edges for column in pair]
    assert all(0 <= column < (1 << rows) for column in columns)
    outside_checks = 0
    for i, (a, b) in enumerate(edges):
        if sperner:
            assert a & ~b, (i, 'first column contained in second')
            assert b & ~a, (i, 'second column contained in first')
        for j, column in enumerate(columns):
            if j // 2 != i:
                assert column & ~(a | b), (i, j, 'covered outside column')
                outside_checks += 1
    return outside_checks


def check_sets(edges):
    """Separate implementation: sets and membership, no bitwise operations."""
    checks = 0
    for i, (a, b) in enumerate(edges):
        assert not a <= b and not b <= a
        for j, pair in enumerate(edges):
            if j != i:
                for column in pair:
                    assert any(x not in a and x not in b for x in column)
                    checks += 1
    return checks


def construction_masks(d):
    assert d >= 3
    r = d // 2
    p, q = 1 << d, 1 << (d + 1)
    edges = [(mask(s) | p, mask(s) | q)
             for s in combinations(range(d), r)]
    u = mask(range(r + 1))
    v = mask(range(r)) | (1 << (r + 1))
    assert u != v and r + 1 < d
    return edges + [(u, v)]


def construction_sets(d):
    r = d // 2
    edges = [(frozenset(s) | {d}, frozenset(s) | {d + 1})
             for s in combinations(range(d), r)]
    return edges + [(frozenset(range(r + 1)),
                     frozenset(range(r)) | {r + 1})]


def main():
    # Verbatim block data from the proof of source Proposition 5.11,
    # displayed by edges rather than by rows of its two-line table.
    source = [({1, 2}, {1, 3}), ({1, 4}, {1, 5}),
              ({2, 4}, {2, 5}), ({3, 4}, {3, 5})]
    source_masks = [(mask(x - 1 for x in a), mask(x - 1 for x in b))
                    for a, b in source]
    fidelity_checks = check_masks(source_masks, 5, True)
    assert check_sets(source) == fidelity_checks == 24
    # Sperner lower bound for eight CFF columns and four ECFF edge unions.
    assert middle(4) < 8 <= middle(5)
    assert middle(3) < 4 <= middle(4)
    ecff4 = [(mask(s), mask(s)) for s in combinations(range(4), 2)][:4]
    assert check_masks(ecff4, 4, False) == 24
    print('FIDELITY Proposition 5.11: t(E_8)=5, t_e(E_8)=4; 24 outside checks')

    results = []
    total = 0
    independent_total = 0
    previous_m = 0
    for d in range(3, 14):
        r = d // 2
        m = middle(d) + 1
        edges = construction_masks(d)
        assert len(edges) == m and m > previous_m
        assert all(a.bit_count() == b.bit_count() == r + 1 for a, b in edges)
        checks = check_masks(edges, d + 2, True)
        assert checks == 2 * m * (m - 1)
        total += checks
        # These are capacity implications of Sperner, not an exhaustive
        # search over all matrices. The universal justification is in README.
        assert middle(d) < m <= middle(d + 1)
        assert middle(d + 1) < 2 * m <= middle(d + 2)
        ecff = [(mask(s), mask(s))
                for s in combinations(range(d + 1), (d + 1) // 2)][:m]
        assert len(ecff) == m
        assert check_masks(ecff, d + 1, False) == checks
        independent = None
        if d <= 7:
            independent = check_sets(construction_sets(d))
            assert independent == checks
            independent_total += independent
        results.append(dict(d=d, r=r, m=m, t_e=d + 1, t=d + 2,
                            outside_checks=checks,
                            independent_set_checks=independent))
        print(f'FAMILY d={d} r={r} m={m} t_e={d+1} t={d+2} outside_checks={checks}')
        previous_m = m

    # Negative control: copying an ordinary block into the exceptional edge
    # violates the outside-column condition and must be rejected.
    broken = construction_masks(3)
    broken[-1] = (broken[0][0], broken[-1][1])
    try:
        check_masks(broken, 5, True)
    except AssertionError:
        print('NEGATIVE_CONTROL rejected a repeated block in distinct edges')
    else:
        raise AssertionError('checker accepted the negative control')
    assert total == 8203100
    print(json.dumps(dict(arxiv='2605.12634v3', fidelity='Proposition 5.11',
                          family=results, total_outside_checks=total,
                          independent_set_checks=independent_total),
                     sort_keys=True, separators=(',', ':')))
    print(f'TOTAL bit-mask CFF outside-column checks={total}; independent set checks={independent_total}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
