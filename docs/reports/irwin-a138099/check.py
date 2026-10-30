#!/usr/bin/env python3
"""Exact finite checks for Irwin's conjecture in OEIS A138099, revision 26.

The hooks are counted from cell coordinates, independently of the triangle and
of the pair-count identity. These bounded checks are not an all-index proof.
"""
from collections import Counter, defaultdict
from math import isqrt

MAX_N = 1001
# OEIS A138099, revision 26, %S/%T/%U; sequence by Paul Curtz.
SOURCE_TERMS = [
    1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7,8,5,6,7,8,9,6,7,8,9,10,6,7,8,9,
    10,11,7,8,9,10,11,12,7,8,9,10,11,12,13,8,9,10,11,12,13,14,8,9,10,11,
    12,13,14,15,9,10,11,12,13,14,15,16,9,10,11,12,13,
]
# The three lines of Irwin's diagram in the %F block (2016-05-14).
SOURCE_DIAGRAM = [
    '                                    o   oo  ooo',
    'Empty, o, oo, oo, oo, oo  , ooo , ooo, ooo, ooo, ...',
    '              o   oo  ooo   ooo   ooo  ooo  ooo',
]


def triangle_terms(limit):
    result = []
    row = 1
    while len(result) < limit:
        result.extend(j + row // 2 for j in range(1, (row + 1) // 2 + 1))
        row += 1
    return result[:limit]


def spiral_points(limit):
    """Start at the origin; walk W,S,E,N, with leg lengths 1,1,2,2,... ."""
    x = y = 0
    if limit:
        yield x, y
    emitted = 1
    leg = 0
    directions = [(-1, 0), (0, 1), (1, 0), (0, -1)]
    while emitted < limit:
        dx, dy = directions[leg % 4]
        for _ in range(leg // 2 + 1):
            x += dx
            y += dy
            yield x, y
            emitted += 1
            if emitted == limit:
                return
        leg += 1


def hook_sum(cells):
    """Literally count self, cells strictly right, and cells strictly below."""
    rows = defaultdict(list)
    columns = defaultdict(list)
    for x, y in cells:
        rows[y].append(x)
        columns[x].append(y)
    return sum(1 + sum(xx > x for xx in rows[y])
               + sum(yy > y for yy in columns[x]) for x, y in cells)


def pair_count(cells):
    rows = Counter(y for x, y in cells)
    columns = Counter(x for x, y in cells)
    return (len(cells) + sum(r * (r - 1) // 2 for r in rows.values())
            + sum(c * (c - 1) // 2 for c in columns.values()))


def normalize(cells):
    if not cells:
        return frozenset()
    min_x = min(x for x, y in cells)
    min_y = min(y for x, y in cells)
    return frozenset((x - min_x, y - min_y) for x, y in cells)


def source_shapes():
    separators = [-1] + [i for i, c in enumerate(SOURCE_DIAGRAM[1]) if c == ',']
    for left, right in zip(separators, separators[1:]):
        yield {(x, y) for y, line in enumerate(SOURCE_DIAGRAM)
               for x in range(left + 1, min(right, len(line))) if line[x] == 'o'}


def main():
    terms = triangle_terms(MAX_N)
    assert terms[:len(SOURCE_TERMS)] == SOURCE_TERMS
    expected_rows = [[1], [2], [2, 3], [3, 4], [3, 4, 5],
                     [4, 5, 6], [4, 5, 6, 7], [5, 6, 7, 8]]
    assert expected_rows == [list(range(1 + row // 2,
                                      1 + row // 2 + (row + 1) // 2))
                             for row in range(1, 9)]
    # Merca's published %F formula, using integer square roots and integer division.
    for n, value in enumerate(terms, 1):
        assert value == n - (isqrt(4 * n - 3) - 1) ** 2 // 4, n
    print(f'FIDELITY: {len(SOURCE_TERMS)} published terms and 8 example rows matched')
    print(f'FIDELITY: Merca formula matched for 1 <= n <= {MAX_N}')

    points = list(spiral_points(MAX_N))
    assert len(set(points)) == MAX_N
    shapes = list(source_shapes())
    assert len(shapes) == 10
    for n, shape in enumerate(shapes):
        assert len(shape) == n, n
        assert normalize(set(points[:n])) == normalize(shape), n
    print('FIDELITY: all 10 source diagrams (empty through 9 cells) matched by translation')

    cells = set()
    totals = [0]
    for n, point in enumerate(points, 1):
        assert point not in cells
        cells.add(point)
        total = hook_sum(cells)
        assert total == pair_count(cells), n
        assert total - totals[-1] == terms[n - 1], n
        totals.append(total)
    assert totals[:10] == [0, 1, 3, 5, 8, 11, 15, 18, 22, 27]
    for n, value in {16: 64, 25: 125, 100: 1000, 1000: 31532}.items():
        assert totals[n] == value, n
    print(f'IDENTITY: direct hook increments matched for 1 <= n <= {MAX_N}')
    print(f'PAIR_COUNT: direct hooks matched for 0 <= n <= {MAX_N}')
    print('H(T_n), n=0..9:', ','.join(map(str, totals[:10])))
    print('First 25 increments:', ','.join(map(str, terms[:25])))
    print('Selected hook sums:', {n: totals[n] for n in [16, 25, 100, 1000, 1001]})

    # All eight rotations/reflections; actual hook directions stay right/down.
    for n in range(101):
        cells = set(points[:n])
        for swap in [False, True]:
            for sx in [-1, 1]:
                for sy in [-1, 1]:
                    image = {(sx * (y if swap else x), sy * (x if swap else y))
                             for x, y in cells}
                    assert hook_sum(image) == totals[n], (n, swap, sx, sy)
    print('ORIENTATION: 8 symmetries matched for 0 <= n <= 100')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
