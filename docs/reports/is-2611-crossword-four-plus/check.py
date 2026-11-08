#!/usr/bin/env python3
"""Exact checks for Isham--Schnitzer, arXiv:2611.03364v1, Conjecture 1.4.

No external packages. True means white. The upper-bound check enumerates an
explicit relaxation containing every Four+ 15-by-15 grid; it does not assume
connectivity, rotational symmetry, or the proposed construction.
"""
from collections import deque
from fractions import Fraction


def runs(line):
    out = []
    start = None
    for j, white in enumerate(list(line) + [False]):
        if white and start is None:
            start = j
        if not white and start is not None:
            out.append((start, j))
            start = None
    return out


def inspect(grid):
    n = len(grid)
    assert n and all(len(row) == n for row in grid)
    row_runs = [runs(row) for row in grid]
    col_runs = [runs([grid[i][j] for i in range(n)]) for j in range(n)]
    assert all(b - a >= 4 for line in row_runs + col_runs for a, b in line)
    assert all(grid[i][j] == grid[n - 1 - i][n - 1 - j]
               for i in range(n) for j in range(n))
    white = {(i, j) for i in range(n) for j in range(n) if grid[i][j]}
    assert white
    seen = {min(white)}
    queue = deque(seen)
    while queue:
        i, j = queue.popleft()
        for cell in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
            if cell in white and cell not in seen:
                seen.add(cell)
                queue.append(cell)
    assert seen == white
    return sum(map(len, row_runs + col_runs)), list(map(len, row_runs)), list(map(len, col_runs))


def construction(k):
    n = 5 * k
    assert k >= 3
    black = set()
    # Independent seed-and-extension specification, not a clue-count formula.
    for i in range(3, n - 3):
        for j in range(3, n - 3):
            if (j - 2 * i) % 5 != 3:
                continue
            black.add((i, j))
            if i == 3:
                black.update((r, j) for r in range(3))
            if i == n - 4:
                black.update((r, j) for r in range(n - 3, n))
            if j == 3:
                black.update((i, c) for c in range(3))
            if j == n - 4:
                black.update((i, c) for c in range(n - 3, n))
    grid = [[(i, j) not in black for j in range(n)] for i in range(n)]
    # Check the alternative clamping formula used in the general argument.
    clamp = lambda x: max(3, min(x, n - 4))
    assert all((not grid[i][j]) == ((clamp(j) - 2 * clamp(i)) % 5 == 3)
               for i in range(n) for j in range(n))
    return grid


def exhaustive_upper_bound():
    # Smallest parameter of the argument: k=3, n=15.
    k, n = 3, 15
    t_mask = sum(1 << j for j in range(n) if (j + 1) % 5 in (2, 3, 4))
    patterns = []
    for mask in range(1 << n):
        line = [bool(mask & (1 << j)) for j in range(n)]
        rr = runs(line)
        if any(b - a < 4 for a, b in rr):
            continue
        r = len(rr)
        assert r <= k
        e = int(line[1]) + int(line[n - 2])
        assert (mask & t_mask).bit_count() >= 2 * r + e
        black_t = t_mask & ~mask
        if r == k:
            assert black_t == 0
            for j in range(1, k):
                assert not (line[5 * j - 1] and line[5 * j])
        patterns.append((k - r, black_t))
    assert len(patterns) == 296
    # best[t]: least sum of column deficits among all pairs that can cover
    # t maximal T-indexed rows. Taking independent minima is a relaxation.
    best = [None] * (3 * k + 1)
    pairs = 0
    for d1, b1 in patterns:
        for d2, b2 in patterns:
            pairs += 1
            s = d1 + d2
            u = (b1 | b2).bit_count()
            assert u <= 2 * k + 2 * s - 2
            if s == 0:
                assert u == 0
            if s == 1:
                assert u <= k + 2
            for t in range(u + 1):
                best[t] = s if best[t] is None else min(best[t], s)
    assert all(s is not None for s in best)
    half_deficit = [3 * k - t + (k - 1) * best[t] for t in range(3 * k + 1)]
    assert min(half_deficit) == 3 * k - 4
    upper = 2 * n * k - 2 * min(half_deficit)
    assert upper == 80 == inspect(construction(k))[0]
    print(f'EXHAUSTIVE n={n}: {1 << n} masks, {len(patterns)} legal lines, {pairs} ordered pairs')
    print(f'EXHAUSTIVE minimum pair deficits by t=0..{3*k}: {best}')
    print(f'EXHAUSTIVE row-T plus column-S deficit bounds: {half_deficit}')
    print(f'EXHAUSTIVE upper bound={upper}; attained={upper}')


def main():
    for n in range(4, 9):
        assert inspect([[True] * n for _ in range(n)])[0] == 2 * n
    # Source Proposition 1.3 and the OEIS example both give a(9)=26.
    black9 = {(0, 4), (1, 4), (7, 4), (8, 4), (4, 0), (4, 1), (4, 7), (4, 8)}
    assert inspect([[(i, j) not in black9 for j in range(9)] for i in range(9)])[0] == 26
    print('FIDELITY all-white n=4..8: 8,10,12,14,16; displayed n=9 grid: 26')
    known = {15: 80, 20: 144, 25: 228, 30: 332}
    for k in range(3, 31):
        n = 5 * k
        count, rows, cols = inspect(construction(k))
        formula = 10 * k * k - 6 * k + 8
        assert Fraction(2 * (n * n - 3 * n + 20), 5) == formula == count
        assert rows.count(k - 1) == cols.count(k - 1) == 3 * k - 4
        assert all(x in (k - 1, k) for x in rows + cols)
        if n in known:
            assert count == known[n]
            print(f'FIDELITY source known value: n={n}, clues={count}')
        if k in (7, 8, 9, 30):
            print(f'CONSTRUCTION n={n}: clues={count}, symmetric=True, connected=True, minimum-run>=4')
    exhaustive_upper_bound()
    # Independent integer enumeration of the stated four-case summation.
    for k in range(3, 101):
        for t in range(3 * k + 1):
            if t == 0:
                need = 0
            elif t <= k + 2:
                need = 1
            elif t <= 2 * k + 2:
                need = 2
            else:
                need = (t - 2 * k + 3) // 2
            assert 3 * k - t + (k - 1) * need >= 3 * k - 4
    print('SUMMATION exact integer cases checked for k=3..100 and every t=0..3k')
    print('SCOPE finite checks support, but do not prove, the unbounded conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
