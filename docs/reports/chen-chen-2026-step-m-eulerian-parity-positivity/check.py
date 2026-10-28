#!/usr/bin/env python3
"""Exact finite checks for Chen--Chen's step-m Eulerian parity conjecture.

The independent implementation uses only Python integers. SOURCE_TERMS contains
rows 1--10 of OEIS A392664, A392665, A392666, respectively (m = -r).
Finite checks support a proposed induction; they do not prove its unbounded claim.
"""

from functools import cache

SOURCE_TERMS = {
    "3": [
        1,
        1,
        1,
        1,
        -4,
        1,
        1,
        3,
        3,
        1,
        1,
        -14,
        -30,
        -14,
        1,
        1,
        17,
        262,
        262,
        17,
        1,
        1,
        -48,
        -1497,
        -4192,
        -1497,
        -48,
        1,
        1,
        79,
        8157,
        50003,
        50003,
        8157,
        79,
        1,
        1,
        -178,
        -42128,
        -514222,
        -1100066,
        -514222,
        -42128,
        -178,
        1,
        1,
        333,
        214200,
        4829952,
        19299834,
        19299834,
        4829952,
        214200,
        333,
        1
    ],
    "4": [
        1,
        1,
        1,
        1,
        -6,
        1,
        1,
        11,
        11,
        1,
        1,
        -44,
        -154,
        -44,
        1,
        1,
        117,
        1562,
        1562,
        117,
        1,
        1,
        -370,
        -12689,
        -34364,
        -12689,
        -370,
        1,
        1,
        1087,
        95853,
        568339,
        568339,
        95853,
        1087,
        1,
        1,
        -3288,
        -695972,
        -8072936,
        -17050170,
        -8072936,
        -695972,
        -3288,
        1,
        1,
        9833,
        4960580,
        104809652,
        409138334,
        409138334,
        104809652,
        4960580,
        9833,
        1
    ],
    "5": [
        1,
        1,
        1,
        1,
        -8,
        1,
        1,
        23,
        23,
        1,
        1,
        -106,
        -414,
        -106,
        1,
        1,
        405,
        5210,
        5210,
        405,
        1,
        1,
        -1644,
        -54585,
        -145880,
        -54585,
        -1644,
        1,
        1,
        6547,
        530721,
        3079435,
        3079435,
        530721,
        6547,
        1,
        1,
        -26222,
        -4966352,
        -55849394,
        -117018530,
        -55849394,
        -4966352,
        -26222,
        1,
        1,
        104849,
        45588716,
        925915724,
        3563737526,
        3563737526,
        925915724,
        45588716,
        104849,
        1
    ]
}


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def triangle(m, last):
    """One-based rows; zero padding makes the displayed recurrence literal."""
    rows = [[], [0, 1]]
    for n in range(2, last + 1):
        row = [0] * (n + 1)
        row[1] = row[n] = 1
        for k in range(2, n):
            row[k] = (m * (n - k) + 1) * rows[n - 1][k - 1] + (
                m * (k - 1) + 1
            ) * rows[n - 1][k]
        rows.append(row)
    return rows


@cache
def direct(m, n, k):
    """Separate recursive evaluation, used only in the small cross-check."""
    if k == 1 or k == n:
        return 1
    return (m * (n - k) + 1) * direct(m, n - 1, k - 1) + (
        m * (k - 1) + 1
    ) * direct(m, n - 1, k)


def main():
    rows_checked = cells_checked = identities = two_steps = direct_checks = 0
    grid = [(r, 403) for r in range(3, 11)] + [(1000000, 30)]
    for r, last in grid:
        rows = triangle(-r, last)
        if str(r) in SOURCE_TERMS:
            actual = [x for n in range(1, 11) for x in rows[n][1:]]
            require(actual == SOURCE_TERMS[str(r)], ('OEIS fidelity', r))
            require(all(x >= 0 for n in range(2, 401, 2)
                        for x in rows[n][1:]), ('source verification', r))
            print(f'FIDELITY m={-r}: 55 published terms matched; '
                  '200 even rows through n=400 entrywise nonnegative')
        product_sum = 1
        for n in range(1, last + 1):
            row = rows[n][1:]
            rows_checked += 1
            cells_checked += n
            require(row == row[::-1], ('reflection', r, n))
            if n >= 2:
                product_sum *= -r * (n - 2) + 2
            require(sum(row) == product_sum, ('source row-sum product', r, n))
            if n % 2 == 0:
                require(all(x > 0 for x in row), ('even strict positivity', r, n))
                require(sum((-1) ** j * x for j, x in enumerate(row)) == 0,
                        ('source alternating sum', r, n))
                require(rows[n][2] >= 1, ('even second-column lower bound', r, n))
            elif n >= 3:
                require(all(x < 0 for x in row[1:-1]), ('odd interior sign', r, n))
            if n >= 2:
                b = rows[n][2]
                require(r * b == 2 * (1 + (-1) ** n * (r - 1) ** (n - 1))
                        - r * (n - 1), ('second-column identity', r, n))
                identities += 1
            if 2 <= n <= last - 2:
                require(rows[n + 2][2] == (r - 1) ** 2 * rows[n][2]
                        + r * (r - 2) * n - r ** 2 + 2,
                        ('second-column two-step identity', r, n))
                two_steps += 1
            if n <= 12:
                for k in range(1, n + 1):
                    require(rows[n][k] == direct(-r, n, k),
                            ('independent recursion', r, n, k))
                    direct_checks += 1
        print(f'GRID r={r}, n=1..{last}: symmetry, row-sum product, '
              'even positivity, odd interior negativity and column identities OK')
    require(rows_checked == 3254, ('row count', rows_checked))
    require(cells_checked == 651713, ('cell count', cells_checked))
    print(f'COUNTS rows={rows_checked} cells={cells_checked} '
          f'closed_column_identities={identities} two_step_identities={two_steps} '
          f'independent_recursive_values={direct_checks}')
    print('SCOPE: finite exact arithmetic only; no Lean proof and no unbounded verification')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
