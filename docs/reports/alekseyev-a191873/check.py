#!/usr/bin/env python3
"""Exact finite checks for the zero-diagonal conjecture in OEIS A191873.

Only integers and exhaustive finite combinatorics are used. These checks do
not establish the assertion for unbounded n and contain no Lean proof.
"""
from itertools import combinations, permutations


def rows_of(mask, n):
    return tuple((mask >> (n * i)) & ((1 << n) - 1) for i in range(n))


def has_rectangle(rows, n):
    return any((rows[i] & rows[j]).bit_count() >= 2
               for i, j in combinations(range(n), 2))


def has_rectangle_by_entries(rows, n):
    # Independent four-cell formulation, including nonadjacent rows/columns.
    return any(all((rows[i] >> j) & 1 for i in rr for j in cc)
               for rr in combinations(range(n), 2)
               for cc in combinations(range(n), 2))


def weight(rows):
    return sum(row.bit_count() for row in rows)


def diagonal_zero(rows, n):
    return all(((rows[i] >> i) & 1) == 0 for i in range(n))


def has_full_row_or_column(rows, n):
    return (any(row == (1 << n) - 1 for row in rows)
            or any(all((row >> j) & 1 for row in rows) for j in range(n)))


def zero_transversal(rows, n):
    return next((p for p in permutations(range(n))
                 if all(((rows[i] >> p[i]) & 1) == 0 for i in range(n))), None)


def permute_columns(rows, p):
    return tuple(sum(((row >> p[j]) & 1) << j for j in range(len(p)))
                 for row in rows)


def exhaustive(n, expected_z, expected_d):
    z = d = -1
    avoiding = lemma_cases = 0
    avoiding_weights = set()
    for mask in range(1 << (n * n)):
        rows = rows_of(mask, n)
        rect = has_rectangle(rows, n)
        assert rect == has_rectangle_by_entries(rows, n), (n, mask, 'predicate')
        if rect:
            continue
        avoiding += 1
        w = weight(rows)
        avoiding_weights.add(w)
        z = max(z, w)
        if diagonal_zero(rows, n):
            d = max(d, w)
        if not has_full_row_or_column(rows, n):
            lemma_cases += 1
            p = zero_transversal(rows, n)
            assert p is not None, (n, mask, 'zero-transversal')
            b = permute_columns(rows, p)
            assert diagonal_zero(b, n) and weight(b) == w
            assert not has_rectangle(b, n)
    assert (z, d) == (expected_z, expected_d), (n, z, d)
    assert avoiding_weights == set(range(z + 1))
    if n >= 2:
        threshold = next(k for k in range(n * n + 1) if k not in avoiding_weights)
        assert threshold == z + 1
    else:
        threshold = None  # A001197 starts at n=2; do not invent its n=1 term.
    extremizers = 0
    if n >= 3:
        assert z == d
        for mask in range(1 << (n * n)):
            if mask.bit_count() != z:
                continue
            rows = rows_of(mask, n)
            if not has_rectangle(rows, n):
                extremizers += 1
                assert not has_full_row_or_column(rows, n)
                p = zero_transversal(rows, n)
                assert p is not None
                b = permute_columns(rows, p)
                assert diagonal_zero(b, n) and weight(b) == z
                assert not has_rectangle(b, n)
    print(f'n={n}: total={1 << (n*n)}, rectangle_free={avoiding}, '
          f'Z={z}, D={d}, forcing_threshold={threshold}, '
          f'zero_transversal_cases={lemma_cases}, '
          f'extremizers_checked={extremizers}')
    return z, d


def main():
    assert has_rectangle((3, 3), 2)
    nonadjacent = (9, 0, 0, 9)
    assert has_rectangle(nonadjacent, 4)
    assert has_rectangle_by_entries(nonadjacent, 4)
    assert not any(all((nonadjacent[i] >> j) & 1
                       for i in (r, r + 1) for j in (c, c + 1))
                   for r in range(3) for c in range(3))
    print('Predicate controls: full 2x2 and nonadjacent rectangle rejected')
    results = [exhaustive(n, z, d) for n, z, d in
               [(1, 1, 0), (2, 3, 2), (3, 6, 6), (4, 9, 9)]]
    assert results[0][0] != results[0][1] and results[1][0] != results[1][1]
    print('Fidelity: A072567 and A191873 terms n=1..4; A001197 terms n=2..4 match')
    print('Excluded-domain controls: equality fails at n=1 and n=2')
    knuth = tuple(int(s[::-1], 2) for s in (
        '01110000', '10010100', '00010011', '01000101',
        '10100010', '11001000', '00101001', '00001110'))
    assert weight(knuth) == 24 and diagonal_zero(knuth, 8)
    assert not has_rectangle_by_entries(knuth, 8)
    print('Fidelity: Knuth A191873 n=8 displayed matrix has weight 24, zero diagonal, no rectangle')
    for n in range(3, 65):
        rows = tuple((1 << i) | (1 << ((i + 1) % n)) for i in range(n))
        assert weight(rows) == 2 * n and not has_rectangle(rows, n)
        p = tuple((j + 2) % n for j in range(n))
        assert len(set(p)) == n
        b = permute_columns(rows, p)
        assert diagonal_zero(b, n) and weight(b) == 2 * n
        assert not has_rectangle(b, n)
    print('Cyclic family: weight 2n, rectangle-free, explicit zero transversal checked for n=3..64')
    print('Scope: exhaustive identity and zero-transversal lemma checked only for n<=4')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
