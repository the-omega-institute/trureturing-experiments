#!/usr/bin/env python3
"""Exact finite checks for the OEIS A102086 row-polynomial conjecture.

Source: https://oeis.org/A102086/internal, revision 19 (2026-07-08).
Python integers and coefficient tuples only; no numerical approximation.
This program supplies finite evidence, not a proof for all rows.
"""
from functools import lru_cache

PUBLISHED = (
    1, 1, 2, 3, 4, 3, 16, 20, 9, 4, 127, 156, 63, 16, 5,
    1363, 1664, 648, 144, 25, 6, 18628, 22684, 8703, 1840,
    275, 36, 7, 311250, 378572, 144243, 29824, 4200, 468, 49, 8,
    6173791, 7504640, 2849400, 582640, 79775, 8316, 735, 64, 9,
    142190703, 172785512,
)


def add(p, q):
    return tuple((p[i] if i < len(p) else 0) +
                 (q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q))))


def scale(c, p):
    return tuple(c * a for a in p)


def triangle(last):
    """NAME recurrence, implemented without using the polynomial table."""
    rows = []
    for n in range(last + 1):
        row = [sum(rows[n - 1][j] * rows[j][k] for j in range(k, n))
               for k in range(n)]
        rows.append(row + [n + 1])
    return rows


@lru_cache(None)
def polynomial(a, b):
    """Literal FORMULA rules, with the source's exterior zero convention."""
    if b < 0 or b > a:
        return (0,)
    if b == 0:
        return (1,)
    if a == b:
        return (0,) * a + (1,)
    return add(polynomial(a, b - 1), scale(b + 1, polynomial(a - 1, b)))


def explicit_paths(a, b):
    """Enumerate individual terminating paths; no recurrence-value memoization."""
    stack = [(a, b, 1)]
    terminal_weights = [0] * (b + 1)
    path_count = 0
    while stack:
        u, v, weight = stack.pop()
        if v == 0:
            terminal_weights[0] += weight
            path_count += 1
        elif u == v:
            terminal_weights[v] += weight
            path_count += 1
        else:
            assert 0 < v < u
            stack.append((u, v - 1, weight))
            stack.append((u - 1, v, weight * (v + 1)))
    return tuple(terminal_weights), path_count


def main():
    last = 30
    rows = triangle(last)
    flat = tuple(a for row in rows for a in row)
    assert flat[:len(PUBLISHED)] == PUBLISHED
    print(f'FIDELITY: all {len(PUBLISHED)} published terms agree')
    square_examples = [[1], [3, 4], [16, 20, 9],
                       [127, 156, 63, 16], [1363, 1664, 648, 144, 25]]
    for n, expected in enumerate(square_examples):
        actual = [sum(rows[n][j] * rows[j][k] for j in range(k, n + 1))
                  for k in range(n + 1)]
        assert actual == expected
    print('FIDELITY: all five published matrix-square example rows agree')
    for n, row in enumerate(rows):
        assert tuple(row) == polynomial(n + 1, n), n
    entries = sum(len(row) for row in rows)
    print(f'IDENTITY: rows 0..{last}, all {entries} entries agree')
    for a in range(last + 2):
        assert polynomial(a, -1) == (0,)
        assert polynomial(a, a + 1) == (0,)
        assert polynomial(a, 0) == (1,)
        assert polynomial(a, a) == (0,) * a + (1,)
    print('BOUNDARIES: negative/exterior columns, zero column and monomial diagonal agree')
    for n in range(1, last + 1):
        rhs = (0,) * n + (n + 1,)
        for j, c in enumerate(polynomial(n, n - 1)):
            rhs = add(rhs, scale(c, polynomial(j + 1, j)))
        assert rhs == polynomial(n + 1, n), n
    print(f'FIRST_RETURN: exact polynomial identity holds for rows 1..{last}')
    for n in range(8):
        actual, paths = explicit_paths(n + 1, n)
        assert actual == polynomial(n + 1, n)
        print(f'PATH_CHECK: row {n}, {paths} individually enumerated paths')
    for n in range(6):
        print(f'ROW {n}: {rows[n]}')
    print('LIMIT: finite checks do not establish the universal conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
