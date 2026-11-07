#!/usr/bin/env python3
"""Exact, finite checks of OEIS A107667 and Kurkov's weighted-ballot identity.

Source: https://oeis.org/A107667, approved revision 25.
All evidence uses Python integers or fractions.Fraction. No external packages.
"""
from fractions import Fraction
from math import factorial

BOUND = 40
SOURCE_ROWS = (
    (1,),
    (4, 2),
    (45, 9, 3),
    (816, 112, 16, 4),
    (20225, 2200, 225, 25, 5),
    (632700, 58176, 4860, 396, 36, 6),
    (23836540, 1920163, 138817, 9408, 637, 49, 7),
    (1048592640, 75683648, 4886464, 290816, 16576, 960, 64, 8),
)
SOURCE_SQUARE_ROWS = (
    (1,),
    (12, 4),
    (216, 45, 9),
    (5248, 816, 112, 16),
    (160675, 20225, 2200, 225, 25),
    (5931540, 632700, 58176, 4860, 396, 36),
    (256182290, 23836540, 1920163, 138817, 9408, 637, 49),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def native_triangle(bound):
    """The source shifted-square equation, without using the conjecture."""
    t = [[0] * (bound + 1) for _ in range(bound + 1)]
    for a in range(bound + 1):
        t[a][a] = a + 1
        for b in range(a - 1, -1, -1):
            t[a][b] = sum(t[a][j] * t[j][b + 1] for j in range(b + 1, a + 1))
    return t


def p_table(m, max_r):
    """The printed base, interior and diagonal rules, including P(1,1,m)=1."""
    p = [[0] * (max_r + 1) for _ in range(max_r + 1)]
    for r in range(1, max_r + 1):
        p[r][1] = 1
        for k in range(2, r):
            p[r][k] = p[r - 1][k] + (r + m) * p[r][k - 1]
        if r >= 2:
            p[r][r] = (r + m) * p[r][r - 1]
    return p


def path_grid(m, bound):
    """Forward propagation of weighted east/north paths constrained by y <= x."""
    grid = [[0] * (bound + 1) for _ in range(bound + 1)]
    grid[0][0] = 1
    for x in range(bound + 1):
        for y in range(x + 1):
            weight = grid[x][y]
            if x < bound:
                grid[x + 1][y] += weight
            if y < x:
                grid[x][y + 1] += weight * (m + x + 1)
    return grid


def enumerate_paths(m, d):
    """Direct traversal of individual paths, independently of all table rules."""
    def visit(x, y, weight):
        if x == d and y == d:
            return weight
        total = 0
        if x < d:
            total += visit(x + 1, y, weight)
        if y < x:
            total += visit(x, y + 1, weight * (m + x + 1))
        return total
    return visit(0, 0, 1)


def check_diagonalization(t, bound):
    """Source T=Q^{-1} D Q, where Q is the earlier, unrelated matrix named P."""
    q = [[Fraction((a + 1) ** (2 * (a - b)), factorial(a - b))
          if b <= a else Fraction(0) for b in range(bound + 1)]
         for a in range(bound + 1)]
    recovered = [[Fraction(0)] * (bound + 1) for _ in range(bound + 1)]
    for a in range(bound + 1):
        for b in range(a + 1):
            recovered[a][b] = ((a + 1) * q[a][b]
                               - sum(q[a][j] * recovered[j][b]
                                     for j in range(b, a)))
            require(recovered[a][b] == t[a][b], ('diagonalization', a, b))
    return (bound + 1) * (bound + 2) // 2


def main():
    t = native_triangle(BOUND)
    for a, row in enumerate(SOURCE_ROWS):
        require(tuple(t[a][:a + 1]) == row, ('source row', a))
    print('FIDELITY: all 36 entries of the 8 explicitly printed T rows agree')
    for a, row in enumerate(SOURCE_SQUARE_ROWS):
        actual = tuple(sum(t[a][j] * t[j][b] for j in range(b, a + 1))
                       for b in range(a + 1))
        require(actual == row, ('source square row', a))
    print('FIDELITY: all 28 entries of the 7 explicitly printed T^2 rows agree')
    count = check_diagonalization(t, 15)
    print(f'FIDELITY: source diagonalization agrees on {count} entries, 0 <= b <= a <= 15')

    grids = {m: path_grid(m, BOUND - m) for m in range(BOUND + 1)}
    entries = 0
    p_entries = 0
    for m in range(BOUND + 1):
        p = p_table(m, BOUND - m + 1)
        grid = grids[m]
        for x in range(BOUND - m + 1):
            for y in range(x + 1):
                require(p[x + 1][y + 1] == grid[x][y], ('P versus paths', m, x, y))
                p_entries += 1
            require(t[x + m][m] == (x + m + 1) * p[x + 1][x + 1],
                    ('conjectured identity', x, m))
            entries += 1
    print(f'IDENTITY: native T, printed P rules and weighted path DP agree on all {entries} entries, 0 <= b <= a <= {BOUND}')
    print(f'PATH INTERPRETATION: {p_entries} P entries agree with forward path DP')

    first_returns = 0
    for m in range(BOUND + 1):
        for d in range(1, BOUND - m + 1):
            rhs = sum((m + ell + 1) * grids[m + 1][ell - 1][ell - 1]
                      * grids[m + ell][d - ell][d - ell]
                      for ell in range(1, d + 1))
            require(grids[m][d][d] == rhs, ('first-return identity', m, d))
            first_returns += 1
    print(f'FIRST RETURN: {first_returns} exact convolution checks')
    direct = 0
    for m in range(6):
        for d in range(9):
            require(enumerate_paths(m, d) == grids[m][d][d], ('direct paths', m, d))
            direct += 1
    print(f'DIRECT ENUMERATION: {direct} checks, 0 <= m <= 5, 0 <= d <= 8')
    print('SCOPE: finite exact evidence only; no all-index or Lean proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
